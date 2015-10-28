import os
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table, Column

import sys
sys.path.append("../dependencies")
import WeightedKDE
from walker_random_sampling import WalkerRandomSampling

from bangs_utils import prepare_data_saving, prepare_plot_saving, \
    BangsDirectories, set_plot_ticks


class PosteriorPredictiveChecks:

    def chi_square(self, y, E_y, sig_y):
        """ 
        An example of a "discrepancy function", the chi-square

        Parameters
        ----------
        y : float
            "Observed" value

        E_y : float
            "Expected" (theoretical) value

        sig_y: float
            Standard deviation of the "observed" values

        Returns 
        ------
        float    

        Notes
        -----
        Negative values of sig_y produce a mask, i.e. those values are not
        considered in the chi-square computation.
        For further details on "discrepancy functions" see Gelman, Meng & Stern (1996).
        """
    
        if y.ndim == 1:
            loc = np.where(sig_y > 0.)[0]
            return np.sum((y[loc]-E_y[loc])**2/sig_y[loc]**2)
        elif y.ndim == 2:
            my = np.ma.masked_array(y, mask=(sig_y <= 0.))
            mE_y = np.ma.masked_array(E_y, mask=(sig_y <= 0.))
            msig_y = np.ma.masked_array(sig_y, mask=(sig_y <= 0.))
            return np.ma.sum((my-mE_y)**2/msig_y**2, axis=0)


    def load(self, file_name):
        """ 
        Load a file containing the posterior predictive checks

        Parameters
        ----------
        file_name : str
            Name of the file.
        """

        name = os.path.join(BangsDirectories.results_dir,
                BangsDirectories.pybangs_data, file_name)

        my_table = Table.read(name)
    
        self.data = my_table

    def compute(self, observed_catalogue, filters, discrepancy=None, 
            n_replicated=2000, file_name=None):
        """ 
        Compute  posterior predictive checks quantities.

        Parameters
        ----------
        observed_catalogue : `bangs_photometry.ObservedCatalogue`
            Class containing an observed photometric catalogue.

        filters : `bangs_filters.PhotometricFilters`
            Class containing a set of photometric filters.

        discrepancy : function, optional
            The discrepancy function used in the posterior predicitve check.

        n_replicated: int, optional
            The number of replicated data to draw.

        file_name : str, optional
            Name of the output catalogue, wuthout including the direcory tree.
            It will be saved into the RESULTS_DIR/PYBANGS_DATA folder (which
            will be created if not present).
        """

        if discrepancy is None:
            discrepancy = self.chi_square

        # Copy from the catalogue the column containing the object IDs
        objID = Column(data=observed_catalogue.data['ID'], name='ID', dtype=np.int32) 

        n_obj = len(observed_catalogue.data['ID'])
        
        # Defines columns containing the number of photometric bands actually
        # used in the BANGS run, for a given object,
        n_used_bands = Column(name='n_used_bands', dtype=np.int32, length=n_obj)

        deg_of_freedom = Column(name='dof', dtype=np.int32, length=n_obj)

        # the average chi-square, and average chi-square/n_used_bands
        # Compute the "average chi square", a measure of the predicitve accuacy
        # of the model (e.g. see Section 6.5 of "Bayesian Data Analysis", by
        # Gelman, Carlin, Stern and Rubin )
        aver_chi_square = Column(name='aver_chi_square', dtype=np.float32,
                length=n_obj)

        # Int from 0 to x of chi^2(x) with N-1 degreed of freedom (see Johnson,
        # V. E. (2004). A Bayesian chi2 Test for Goodness-of-Fit on JSTOR.
        # Annals of Statistics for an explanation of why the average chi^2 has
        # N-1 and not N-k-1 degrees of freedom)

        left_cumul_probability = Column(name='left_cumul_probability', dtype=np.float32,
                length=n_obj)

        # Int from x to +infty of chi^2(x)
        right_cumul_probability = Column(name='right_cumul_probability', dtype=np.float32,
                length=n_obj)

        aver_red_chi_square = Column(name='aver_red_chi_square',
                dtype=np.float32, length=n_obj)

        p_value = Column(name='p_value',
                dtype=np.float32, length=n_obj)

        my_cols = [objID, n_used_bands, deg_of_freedom, aver_chi_square,
                aver_red_chi_square, left_cumul_probability,
                right_cumul_probability, p_value]

        my_table = Table(my_cols)

        obs_flux = np.zeros(filters.n_bands, np.float32)
        obs_flux_err = np.zeros(filters.n_bands, np.float32)

        model_flux = np.zeros(filters.n_bands, np.float32)

        jy = 1.E-26

        for i in range(n_obj):

            strID = str(objID[i])
            file = os.path.join(BangsDirectories.results_dir, strID + "_BANGS.fits.gz")

            if os.path.isfile(file):

                # Open the FITS file containing BANGS results for the current object
                hdulist = fits.open(file)
                bangs_data = hdulist['MARGINAL PHOTOMETRY'].data

                # Load the posterior probability and create array of row indices
                probability = hdulist['POSTERIOR PDF'].data['probability']
                n_samples = len(probability)
                row_indices = np.arange(n_samples)

                # Now, draw the weighted samples with replacement
                wrand = WalkerRandomSampling(probability, keys=row_indices)
                replic_data_rows = wrand.random(n_replicated)

                obs_flux = np.zeros(filters.n_bands, np.float32)
                obs_flux_err = np.zeros(filters.n_bands, np.float32)
                model_flux = np.zeros((filters.n_bands, n_samples), np.float32)

                # The replicated data are just the fluxes predicted by your
                # model, drawn from the posterior probability distribution
                # accordingly to their probability, with the effect of
                # observational noise added
                noiseless_flux = np.zeros((filters.n_bands, n_replicated), np.float32)
                replic_flux = np.zeros((filters.n_bands, n_replicated), np.float32)

                n_data = 0

                for j in range(filters.n_bands):

                    # observed flux and its error
                    name = filters.data['flux_colName'][j]
                    obs_flux[j] = observed_catalogue.data[i][name] * filters.units / jy

                    name = filters.data['flux_errcolName'][j]
                    obs_flux_err[j] = observed_catalogue.data[i][name] * filters.units  / jy
                    has_measure = False
                    if obs_flux_err[j] > 0.:
                        has_measure = True

                    # model flux
                    name = '_' + filters.data['label'][j] + '_'
                    model_flux[j,:] = bangs_data[name] / jy

                    if has_measure > 0.:
                        # if defined, add the minimum error in quadrature
                        obs_flux_err[j] = (np.sqrt((obs_flux_err[j]/obs_flux[j])**2 +
                            np.float32(filters.data['min_rel_err'][j])**2) *
                            obs_flux[j])
                        n_data += 1
                        
                # You save in this array the noise-less flux predicted by the model        
                noiseless_flux = model_flux[:, replic_data_rows]

                for j in range(filters.n_bands):
                    # Here you add the observational noise to the
                    # noise-less fluxes predicted by the model, obtaining
                    # the actual replicated data
                    if obs_flux_err[j] > 0.:
                        replic_flux[j,:] = noiseless_flux[j, :] +   \
                        np.random.normal(scale=obs_flux_err[j], size=n_replicated)
                    else:    
                        replic_flux[j,:] = -99.99999

                # Write the replicated data to an output FITS file
                # Create the new FITS file
                new_hdu = fits.HDUList(fits.PrimaryHDU())

                # Add column allowing you to match each row of the replicated
                # data to the row of noise-less fluxes predicted by your model,
                # i.e. those in the "MARGINAL PHOTOMETRY" extension of the
                # BANGS output FITS file 
                # NB: the column indexing start with 1 !!
                ID_col = fits.Column(name='row_index', format='I')

                # Copy the columns defined in the "MARGINAL PHOTOMETRY"
                # extension to the new FITS file
                cols = hdulist['MARGINAL PHOTOMETRY'].columns

                new_hdu.append(fits.BinTableHDU.from_columns(ID_col + cols, nrows=n_replicated, fill=True))

                # Fill with the replicated data fluxes
                j = 0
                for col_name in new_hdu[1].columns.names:
                    if col_name in cols.names:
                        new_hdu[1].data[col_name] = replic_flux[j,:]
                        j += 1

                # Fill with the row indices
                new_hdu[1].data['row_index'] = replic_data_rows

                # Save the file
                name = strID + "_BANGS_replic_data.fits.gz"
                name = prepare_data_saving(name)
                new_hdu.writeto(name)

                # Extend the arrays containing the observed flux and its
                # error to match the shape of the replicated data array
                ext_obs_flux = obs_flux.reshape(filters.n_bands, 1).repeat(n_replicated, 1)
                ext_obs_flux_err = obs_flux_err.reshape(filters.n_bands, 1).repeat(n_replicated, 1)

                # Compute the "discrepancy" for the actual data, and for the replicated data
                discrepancy_data = discrepancy(ext_obs_flux, noiseless_flux, ext_obs_flux_err)
                discrepancy_repl_data = discrepancy(replic_flux, noiseless_flux, ext_obs_flux_err)

                # The p-value is just the fraction of objects for which the
                # discrepancy for the replicated data is larger then that for
                # the actual data! 
                my_table['p_value'][i] = 1. * np.count_nonzero((discrepancy_repl_data > discrepancy_data)) / n_replicated

                dof = n_data
                my_table['n_used_bands'][i] = n_data
                my_table['dof'][i] = dof

                # Here you cosider all samples in the posterior
                ext_obs_flux = obs_flux.reshape(filters.n_bands, 1).repeat(n_samples, 1)
                ext_obs_flux_err = obs_flux_err.reshape(filters.n_bands, 1).repeat(n_samples, 1)
                av_chi_square = np.sum(probability*self.chi_square(ext_obs_flux, model_flux, ext_obs_flux_err)) / np.sum(probability)
                my_table['aver_chi_square'][i] = av_chi_square
                my_table['aver_red_chi_square'][i] = av_chi_square / dof

                cdf = stats.chi2.cdf(av_chi_square, dof)
                my_table['left_cumul_probability'][i] = cdf
                my_table['right_cumul_probability'][i] = 1.-cdf

                hdulist.close()

        self.columns = my_cols
        self.data = my_table

        if file_name is not None:
            name = prepare_data_saving(file_name)
            my_table.write(name)

    def plot_chi2(self, 
            plot_name="BANGS_average_chi_square.pdf"):
        """ 
        Plots the distribution (histogram) of the average chi-square.

        Parameters
        ----------
        plot_name: str, optional
            File name of the output plot.
        """ 

        xdata = self.data['aver_chi_square']
        n_data = len(xdata)
        min_x = 0.
        max_x = 50.

        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)

        ax.set_ylabel("Number of galaxies")
        ax.set_xlabel("$\langle \chi^2 \\rangle$")

        ax.set_xlim((min_x, max_x))

        # Set the correct number of major and mnor tick marks
        set_plot_ticks(ax, prune_y='lower')

        # Plot the histogram of the average chi-square
        kwargs = {'alpha':0.7, 'linewidth':0.5}
        n, bins, patches = ax.hist(xdata, 
                bins=50, 
                range=(min_x, max_x),
                color='gray',
                **kwargs)

        # Integral
        dx = bins[1:]-bins[0:-1]
        norm = np.sum(n*dx)

        # Now, compute the theoretical distribution, i.e. the sum of different
        # chi-square distribution for the different degrees of freedom
        mask = np.zeros(n_data, dtype=bool)
        mask[self.data['dof'] > 0] = True

        min_dof = np.amin(self.data['dof'][mask])
        max_dof = np.amax(self.data['dof'][mask])
                
        dof_range = np.arange(min_dof, max_dof+1)
        frac_data = np.zeros(len(dof_range), dtype=np.float32)

        for i in range(len(dof_range)):
            dof = dof_range[i]
            loc = np.where(self.data['dof'] == dof)[0]
            frac_data[i] = 1.*len(loc)/n_data

        frac_data /= np.sum(frac_data)

        xdata = np.linspace(min_x, max_x, 1000)
        chi_distr = np.zeros(len(xdata), dtype=np.float32)
        for i in range(len(dof_range)):
            dof = dof_range[i]
            chi_distr += frac_data[i] * stats.chi2.pdf(xdata, dof)

        chi_distr *= norm
        ax.plot(xdata, 
                chi_distr,
                color='black',
                linestyle='--')

        name = prepare_plot_saving(plot_name)
        fig.savefig(name, dpi=None, facecolor='w', edgecolor='w',
                orientation='portrait', papertype='a4', format="pdf",
                transparent=False, bbox_inches="tight", pad_inches=0.1)

        plt.close(fig)

    def plot_p_value(self, 
            plot_name="BANGS_p_value.pdf"):
        """ 
        Plots the distribution (histogram) of the p-value.

        Parameters
        ----------
        plot_name: str, optional
            File name of the output plot.
        """ 

        xdata = np.log10(self.data['p_value'])
        n_data = len(xdata)
        min_x = 0.001
        max_x = 0.

        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)

        ax.set_ylabel("Number of galaxies")
        ax.set_xlabel("$\log p$-value")

        ax.set_xlim((min_x, max_x))

        # Set the correct number of major and mnor tick marks
        set_plot_ticks(ax, n_x=3, prune_y='lower')

        # Plot the histogram of the average chi-square
        kwargs = {'alpha':0.7, 'linewidth':0.5}
        n, bins, patches = ax.hist(xdata, 
                bins=50, 
                range=(min_x, max_x),
                color='gray',
                **kwargs)

        y0, y1 = ax.get_ylim()
        levels = (0.01,)
        for lev in levels:
            l = np.log10(lev)
            frac = 1.*len(np.where(xdata <= l)[0])/n_data
            print "Fraction of galaxies with p-value < " + "{:.2f}".format(lev) + " = {:.2f}".format(frac)
            ax.plot((l, l),
                    (y0, y1),
                    color='black',
                    linestyle='--')

        name = prepare_plot_saving(plot_name)
        fig.savefig(name, dpi=None, facecolor='w', edgecolor='w',
                orientation='portrait', papertype='a4', format="pdf",
                transparent=False, bbox_inches="tight", pad_inches=0.1)

        plt.close(fig)

