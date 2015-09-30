import os
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table, Column

import sys
sys.path.append("../dependencies")
import WeightedKDE

from bangs_utils import prepare_data_saving, prepare_plot_saving, \
    BangsDirectories, set_plot_ticks

class PosteriorPredictiveChecks:


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

    def compute(self, observed_catalogue, filters, 
            file_name=None):
        """ 
        Compute  posterior predictive checks quantities.

        Parameters
        ----------
        observed_catalogue : `bangs_photometry.ObservedCatalogue`
            Class containing an observed photometric catalogue.

        filters : `bangs_filters.PhotometricFilters`
            Class containing a set of photometric filters.

        file_name : str, optional
            Name of the output catalogue, wuthout including the direcory tree.
            It will be saved into the RESULTS_DIR/PYBANGS_DATA folder (which
            will be created if not present).
        """

##        if filters is not None:
##            self.filters = filters
##
##        if not hasattr(self, 'filters'):
##            raise AttributeError("No 'PhotometricFilters' class has been passed"
##                    " to the class constructor or to this function!")
##
##        if results_dir is not None:
##            self.results_dir = results_dir
##
##        if not hasattr(self, 'results_dir'):
##            raise AttributeError("No 'results_dir' defined!")

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

        my_cols = [objID, n_used_bands, deg_of_freedom, aver_chi_square,
                aver_red_chi_square, left_cumul_probability,
                right_cumul_probability]

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

                probability = hdulist['POSTERIOR PDF'].data['probability']

                n_samples = len(bangs_data.field(0))
                chi_square = np.zeros(n_samples, np.float32)
                n_data = 0

                for j in range(filters.n_bands):

                    # observed flux and its error
                    name = filters.data['flux_colName'][j]
                    obs_flux = observed_catalogue.data[i][name] * filters.units / jy

                    name = filters.data['flux_errcolName'][j]
                    obs_flux_err = observed_catalogue.data[i][name] * filters.units  / jy
                    has_measure = False
                    if obs_flux_err > 0.:
                        has_measure = True

                    # model flux and its error
                    name = '_' + filters.data['label'][j] + '_'
                    model_flux = bangs_data[name] / jy

                    if has_measure > 0.:
                        # if defined, add the minimum error in quadrature
                        obs_flux_err = (np.sqrt((obs_flux_err/obs_flux)**2 +
                            np.float32(filters.data['min_rel_err'][j])**2) *
                            obs_flux)
                        n_data += 1
                        chi_square += ((obs_flux-model_flux) / obs_flux_err)**2

                dof = n_data

                my_table['n_used_bands'][i] = n_data
                my_table['dof'][i] = dof

                av_chi_square = np.sum(probability*chi_square) / np.sum(probability)
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

        xdata = np.log10(self.data['right_cumul_probability'])
        n_data = len(xdata)
        min_x = -3.
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

