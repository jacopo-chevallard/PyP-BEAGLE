from __future__ import absolute_import
from __future__ import print_function
import logging
import os
from scipy.interpolate import interp1d
from collections import OrderedDict
from bisect import bisect_left
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as plticker
#import pandas as pd
# SEABORN creates by default plots with a filled background!!
#import seaborn as sns
from astropy.io import ascii
from astropy.io import fits

import sys
import pyp_beagle.dependencies.WeightedKDE 
from pyp_beagle.dependencies.walker_random_sampling import WalkerRandomSampling

from .beagle_utils import BeagleDirectories, prepare_plot_saving, set_plot_ticks, \
        prepare_violin_plot, plot_exists, pause, extract_row, is_FITS_file
from .beagle_filters import PhotometricFilters
from .beagle_summary_catalogue import BeagleSummaryCatalogue
from .beagle_residual_photometry import ResidualPhotometry
from .beagle_multinest_catalogue import MultiNestCatalogue
from .beagle_posterior_predictive_checks import PosteriorPredictiveChecks
from .beagle_observed_catalogue import ObservedCatalogue
from six.moves import range


Jy = np.float32(1.E-23)
microJy = np.float32(1.E-23 * 1.E-06)
nanoJy = np.float32(1.E-23 * 1.E-09)
c_light = 2.99792e+18 # Ang/s

p_value_lim = 0.05



class PhotometricCatalogue(ObservedCatalogue):

    def extract_fluxes(self, filters, ID, key='ID', aper_corr=1.):
        """ 
        Extract fluxes and error fluxes for a single object (units are Jy).

        Parameters
        ----------
        filters : class
            Contains the photometric filters

        ID : int, str
            Contains the object ID

        Returns    
        -------
        flux : array
            In units of Jy

        flux_error : array 
            In units of Jy

        Notes
        -----
        The routine also adds in quadrature the minimum relative error defined int he filters class.

        """

        flux = np.zeros(filters.n_bands, np.float32)
        flux_err = np.zeros(filters.n_bands, np.float32)

        row = extract_row(self.data, ID, key=key)

        for j in range(filters.n_bands):

            # observed flux and its error
            name = filters.data['flux_colName'][j]
            if not name:
                flux[j] = -99.
                flux_err[j] = -99.
                continue

            flux[j] = row[name] * aper_corr * filters.units / Jy

            name = filters.data['flux_errcolName'][j]
            flux_err[j] = row[name] * aper_corr * filters.units / Jy

            if flux_err[j] > 0.:
                # if defined, add the minimum error in quadrature
                flux_err[j] = (np.sqrt((flux_err[j]/flux[j])**2 +
                    np.float32(filters.data['min_rel_err'][j])**2) *
                    abs(flux[j]))

        return flux, flux_err

class Photometry:

    def __init__(self, **kwargs):
        
        self.inset_fontsize = BeagleDirectories.inset_fontsize_fraction * BeagleDirectories.fontsize

        self.filters = PhotometricFilters()

        self.observed_catalogue = PhotometricCatalogue()

        self.summary_catalogue = BeagleSummaryCatalogue()

        self.multinest_catalogue = MultiNestCatalogue()

        self.residual = ResidualPhotometry()

        self.PPC = PosteriorPredictiveChecks()

        self.key = kwargs.get('ID_key', 'ID')
        
        self.x_log = kwargs.get('plot_log_wl', False)

        self.plot_full_SED = kwargs.get('plot_full_SED', False)

        self.plot_MAP_SED = kwargs.get('plot_MAP_SED', False)

        self.log_flux = kwargs.get('log_flux', False)

        self.plot_filter_labels = kwargs.get('plot_filter_labels', False)

        self.flux_units = kwargs.get('flux_units', 'microJy')

        self.single_solutions = None
        if kwargs.get('plot_single_solution') is not None:
            self.single_solutions = OrderedDict()
            with fits.open(kwargs.get('plot_single_solution')) as f:
                self.single_solutions['ID'] = f[1].data['ID']
                self.single_solutions['row'] = f[1].data['row_index']

    def plot_marginal(self, ID, max_interval=99.7, 
            print_text=False, print_title=False, replot=False, show=False, units='nanoJy',
            SED_prob_log_scale=False, n_SED_to_plot=10):
        """ 
        Plot the fluxes predicted by BEAGLE.

        The fluxes here considered are those predicted by BEAGLE, given the
        posterior distribution of the model parameters. These are *not*
        replicated data.

        Parameters
        ----------
        ID : int
            ID of the galaxy whose marginal photometry will be plotted.

        max_interval : float, optional
            The marginal photometry is shown to include `max_interval`
            probability, e.g. `max_interval` = 68. will show the 68 % (i.e.
            '1-sigma') (central) credible region of the marginal photometry.

        print_text : bool, optional
            Whether to print further information on the plot, such as
            chi-square, p-value, or leave it empty and neat.

        print_text : bool, optional
            Whether to print the object ID on the top of the plot.

        replot: bool, optional
            Whether to redo the plot, even if it already exists
        """

        if self.flux_units == 'milliJy':
            flux_factor = 1.E+03
            ylabel = "$f_{\\nu}/\\textnormal{mJy}$"
        elif self.flux_units == 'microJy':
            flux_factor = 1.E+06
            ylabel = "$f_{\\nu}/\\upmu\\textnormal{Jy}$"
        elif self.flux_units == 'nanoJy':
            flux_factor = 1.E+09
            ylabel = "$f_{\\nu}/\\textnormal{nJy}$"
        else:
            raise ValueError("Flux units `" + self.flux_units + "` not recognised!")

        # Name of the output plot
        plot_name = str(ID)+'_BEAGLE_marginal_SED_phot.pdf'

        # Check if the plot already exists
        if plot_exists(plot_name) and not replot and not show:
            logging.warning('The plot "' + plot_name + '" already exists. \n Exiting the function.')
            return

        # From the (previously loaded) observed catalogue select the row
        # corresponding to the input ID
        observation = extract_row(self.observed_catalogue.data, ID, key=self.key)

        # Check if you need to apply an aperture correction to the catalogue fluxes
        if 'aper_corr' in self.observed_catalogue.data.dtype.names:
            aper_corr = 10.**(-0.4*observation[0]['aper_corr'])
        else:
            aper_corr = 1.

        # Put observed photometry and its error in arrays
        obs_flux, obs_flux_err = self.observed_catalogue.extract_fluxes(self.filters, ID, key=self.key)
        obs_flux *= flux_factor
        obs_flux_err *= flux_factor

        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)

        # Open the file containing BEAGLE results
        fits_file = os.path.join(BeagleDirectories.results_dir,
                str(ID) + '_' + BeagleDirectories.suffix + '.fits.gz')

        hdulist = fits.open(fits_file)

        # Consider only the extension containing the predicted model fluxes
        old_API = False
        try:
            model_sed = hdulist['marginal photometry']
            old_API = True
        except:
            model_sed = hdulist['apparent magnitudes']

        probability = hdulist['posterior pdf'].data['probability']

        n_bands = len(obs_flux)
        median_flux = np.zeros(n_bands)
        pdf_norm = np.zeros(n_bands)
        _max_y = np.zeros(n_bands)
        min_flux = np.zeros(n_bands)
        max_flux = np.zeros(n_bands)

        y_plot = list(range(n_bands))
        x_plot = list(range(n_bands))

        has_pdf = list(range(n_bands))
        kde_pdf = list(range(n_bands))
        nXgrid = 1000

        if self.x_log:
            wl_eff = np.log10(self.filters.data['wl_eff'])
        else:
            wl_eff = np.array(self.filters.data['wl_eff'])

        # Sort wl_eff array
        sor = np.argsort(wl_eff)

        obs_flux, obs_flux_err, wl_eff = obs_flux[sor], obs_flux_err[sor], wl_eff[sor]

        ok = np.where(obs_flux_err > 0.)[0]

        width = 5*np.min(wl_eff[1:]-wl_eff[0:-1])

        kwargs = {'color':'tomato', 'alpha':0.8, 'edgecolor':'black', 'linewidth':0.2}

        for i in range(n_bands):

            if old_API:
                band_name = self.filters.data['label'][sor[i]]
                xdata = model_sed.data['_'+band_name+'_'] / Jy * flux_factor
            else:
                band_name = self.filters.data['name'][sor[i]] + "_APP"
                xdata = 10.**(0.4*(8.9-model_sed.data[band_name])) * flux_factor

            min_x = np.min(xdata)
            max_x = np.max(xdata)

            min_flux[i] = np.min([min_x, obs_flux[i]-obs_flux_err[i]])
            max_flux[i] = np.max([max_x, obs_flux[i]+obs_flux_err[i]])

            # if min_x == max_x, then you can not use weighted KDE, since you
            # just have one value for the x...this usually happens bacause of
            # IGM absorption, which absorbs the flux blue-ward 1216 AA, making
            # all flux = 0
            if min_x == max_x:
                has_pdf[i] = False
                median_flux[i] = min_x
                continue

            # Compute the marginal PDF through a weighted KDE
            has_pdf[i] = True

            # This function provides you with all the necessary info to draw violin plots
            kde_pdf[i], pdf_norm[i], median_flux[i], x_plot[i], y_plot[i] = prepare_violin_plot(xdata, weights=probability) 

            _max_y[i] = np.max(y_plot[i])

        delta_wl = wl_eff[1:]-wl_eff[0:-1]
        delta_wl = np.concatenate(([delta_wl[0]], delta_wl))
        delta_wl /=  2.

        for i in range(n_bands):

            dwl = delta_wl[i]
            if i > 1:
                dwl = np.min(delta_wl[i-1:i])

            if has_pdf[i]:

                w = 0.4 * dwl / _max_y[i]

                y_grid = np.full(len(x_plot[i]), wl_eff[i])
                _lim_y = kde_pdf[i](median_flux[i])/pdf_norm[i] * w

                ax.fill_betweenx(x_plot[i],
                        y_grid - y_plot[i]*w,
                        y_grid + y_plot[i]*w,
                        **kwargs
                        )

                ax.plot([wl_eff[i]-_lim_y, wl_eff[i]+_lim_y],
                        [median_flux[i], median_flux[i]],
                        color = 'black',
                        linewidth = 0.2
                        )

            ax.plot(wl_eff[i],
                    median_flux[i],
                    color = 'black',
                    marker = "o",
                    markersize = 5,
                    zorder=3,
                    alpha = 0.6
                    )


        # Plot the full SED
        if 'full sed wl' in hdulist: 

            if self.plot_full_SED or self.plot_MAP_SED:

                _n_SED_to_plot = 0
                if self.plot_full_SED: 
                    _n_SED_to_plot = n_SED_to_plot

                if self.plot_MAP_SED:
                    _n_SED_to_plot += 1

                wl = hdulist['full sed wl'].data['wl'][0,:]
                redshifts = hdulist['galaxy properties'].data['redshift']

                if SED_prob_log_scale:
                    max_prob = np.log10(np.amax(probability))
                    min_prob = np.log10(np.amin(probability))
                else:
                    max_prob = np.amax(probability)
                    min_prob = np.amin(probability)

                indices = np.arange(len(probability))

                wrand = WalkerRandomSampling(probability, keys=indices)
                rand_indices = wrand.random(n_SED_to_plot)

                for j in range(_n_SED_to_plot):

                    if self.plot_MAP_SED and j == 0:
                        i = np.argmax(probability)
                        color = "black"
                        lw=1.2
                        max_alpha = 0.8
                    else:
                        i = rand_indices[j]
                        color = "black"
                        lw=0.5
                        max_alpha = 0.4

                    SED = hdulist['full sed'].data[i,:]

                    z = 0.
                    if redshifts[i] > 0.:
                        z = redshifts[i]

                    # Redshift the SED and wl
                    flux_obs = SED / (1.+z)
                    wl_obs = wl * (1.+z)

                    # Convert F_lambda [erg s^-1 cm^-2 A^-1] ----> F_nu [erg s^-1 cm^-2 Hz^-1]
                    flux_obs = (wl_obs)**2/c_light*flux_obs

                    # Scale to nanoJy
                    flux_obs = flux_obs / Jy * flux_factor

                    prob = probability[i]
                    if SED_prob_log_scale:
                        alpha = (np.log10(prob)-min_prob)/(max_prob-min_prob)
                    else:
                        alpha = (prob-min_prob)/(max_prob-min_prob)

                    if self.x_log:
                        wl_obs = np.log10(wl_obs)

                    alpha=1.
                    ax.plot(wl_obs, 
                            flux_obs,
                            color=color,
                            ls="-",
                            lw=lw,
                            alpha=alpha*max_alpha)


        # Determine min and max values of y-axis
        yMax = np.max(max_flux)
        yMin = np.min(np.concatenate((obs_flux[ok], min_flux)))

        dY = yMax-yMin

        yMax += dY * 0.1
        if self.log_flux:
            yMin = 0.
        else:
            yMin -= dY * 0.1

        ax.set_ylim([yMin, yMax])

        x0 = wl_eff[0]
        x1 = wl_eff[-1]
        dx = x1-x0
        ax.set_xlim([x0-0.05*dx, x1+0.05*dx])

        x0, x1 = ax.get_xlim()
        if yMin < 0.: plt.plot( [x0,x1], [0.,0.], color='gray', lw=0.8 )

        # Plot labels of photometric filters
        if self.plot_filter_labels:
            for i in range(n_bands):
                x = wl_eff[i]
                y0, y1 = ax.get_ylim() ; dy = y1-y0
                y = np.min([max_flux[i]+dy*0.2, y1-0.05*dy])

                if old_API:
                    band_name = self.filters.data['label'][sor[i]]
                else:
                    band_name = self.filters.data['name'][sor[i]]

                band_name.replace("_", "") ; band_name = "$" + band_name + "$"

                rotation = 0
                if len(band_name) >= 5:
                    rotation = 45

                ax.text(x, y, 
                        band_name, 
                        fontsize=self.inset_fontsize, 
                        rotation=rotation,
                        color='black',
                        va='center',
                        ha='center')

                ydots = y-0.025*dy
                y = max_flux[i]
                ax.plot([x,x], [y,ydots], 
                        ls=":",
                        color='darkgrey',
                        lw=1.0,
                        zorder=0)


        # Define plotting styles
        if self.x_log:
            ax.set_xlabel("$\\log (\lambda_\\textnormal{eff} / \\textnormal{\AA})$")
        else:
            ax.set_xlabel("$\\lambda_\\textnormal{eff} / \\textnormal{\AA}$")

        ax.set_ylabel(ylabel)

        # Set better location of tick marks
        set_plot_ticks(ax, n_x=5)

        kwargs = {'alpha':0.7}

        plt.errorbar(wl_eff[ok], 
                obs_flux[ok], 
                yerr = obs_flux_err[ok],
                color = "dodgerblue",
                ls = " ",
                marker = "D",
                markeredgewidth = 0.,
                markersize = 8,
                elinewidth=1.0,
                capsize=3,
                **kwargs)

        if self.single_solutions is not None:
            row =  self.single_solutions['row'][self.single_solutions['ID']==ID]
            solution = np.zeros(n_bands, dtype=np.float32)
            for i, band_name in enumerate((self.filters.data['label'][sor])):
                solution[i] = model_sed.data['_'+band_name+'_'][row] / nanoJy

            ax.plot(wl_eff,
                    solution,
                    color = 'green',
                    marker = "*",
                    ls="",
                    markersize = 10,
                    alpha = 0.7
                    )

        if self.log_flux:
            ax.set_yscale('symlog')
            which = 'x'

        # Title of the plot is the object ID
        if print_title: plt.title(str(ID))

        # Location of printed text
        x0, x1 = ax.get_xlim()
        x = x0 + (x1-x0)*0.03
        y0, y1 = ax.get_ylim()
        y = y1 - (y1-y0)*0.10

        if print_text:

            # Print the evidence
            try:
                ax.text(x, y, "$\\log(Z)=" + "{:.2f}".format(self.logEvidence) + "$", fontsize=10 )
            except AttributeError:
                print("ciao")

            # Print the average reduced chi-square
            try:
                row = extract_row(self.PPC.data, ID, key=self.key)
                aver_chi_square = row['aver_chi_square']
                y = y1 - (y1-y0)*0.15
                ax.text(x, y, "$\\langle\chi^2\\rangle=" + "{:.2f}".format(aver_chi_square) + "$", fontsize=10 )
            except AttributeError:
                print("`PosteriorPredictiveChecks` not computed/loaded, hence " \
                "<chi^2> for the object `" + str(ID) + "` is not available")

            try:
                row = extract_row(self.PPC.data, ID, key=self.key)
                aver_red_chi_square = row['aver_red_chi_square']
                n_data = row['n_used_bands']
                y = y1 - (y1-y0)*0.20
                ax.text(x, y,
                        "$\\langle\chi^2/(\\textnormal{N}_\\textnormal{data}-1)\\rangle=" \
                        + "{:.2f}".format(aver_red_chi_square) + "\; \
                        (\\textnormal{N}_\\textnormal{data}=" + \
                        "{:d}".format(n_data) + ")" + "$", fontsize=10 )
            except AttributeError:
                print("`PosteriorPredictiveChecks` not computed/loaded, hence " \
                "<chi^2_red> for the object `" + str(ID) + "` is not available")

        if y0 < 0.: plt.plot( [x0,x1], [0.,0.], color='gray', lw=1.0 )

        if show:
            plt.show()
        else:
            name = prepare_plot_saving(plot_name)

            fig.savefig(name, dpi=None, facecolor='w', edgecolor='w',
                    orientation='portrait', papertype='a4', format="pdf",
                    transparent=False, bbox_inches="tight", pad_inches=0.1)

        plt.close(fig)

        hdulist.close()

    def plot_replicated_data(self, ID, max_interval=99.7, n_replic_to_plot=16,
            print_text=False, replot=False):    
        """ 
        Plot the replicated data.

        Parameters
        ----------
        ID : int
            ID of the galaxy whose marginal photometry will be plotted.

        max_interval : float, optional
            The marginal photometry is shown to include `max_interval`
            probability, e.g. `max_interval` = 68. will show the 68 % (i.e.
            '1-sigma') (central) credible region of the marginal photometry.

        n_replic_to_plot: int, optional
            The number of replicated data that will be plotted. It can be given
            as a single number, or as a pair (n_x, n_y), in which case the
            total number of replicated data plotted will be n_x * n_y

        print_text : bool, optional
            Whether to print further information on the plot, such as
            chi-square, p-value, or leave it empty and neat.

        replot: bool, optional
            Whether to redo the plot, even if it already exists
        """

        if self.flux_units == 'milliJy':
            flux_factor = 1.E+03
            ylabel = "$f_{\\nu}/\\textnormal{mJy}$"
        elif self.flux_units == 'microJy':
            flux_factor = 1.E+06
            ylabel = "$f_{\\nu}/\\upmu\\textnormal{Jy}$"
        elif self.flux_units == 'nanoJy':
            flux_factor = 1.E+09
            ylabel = "$f_{\\nu}/\\textnormal{nJy}$"
        else:
            raise ValueError("Flux units `" + self.flux_units + "` not recognised!")

        # Name of the output plot
        plot_name = str(ID)+'_BEAGLE_replic_data_phot.pdf'

        # Check if the plot already exists
        if plot_exists(plot_name) and not replot:
            logging.warning('The plot "' + plot_name + '" already exists. \n Exiting the function.')
            return

        n_replic_to_plot = np.array(n_replic_to_plot)
        if n_replic_to_plot.size == 1:
            n_plot_x = int(np.sqrt(n_replic_to_plot))
            n_plot_y = n_plot_x
        else:
            n_plot_x = n_replic_to_plot[0]
            n_plot_y = n_replic_to_plot[1]

        # From the (previously loaded) observed catalogue select the row
        # corresponding to the input ID
        observation = extract_row(self.observed_catalogue.data, ID, key=self.key)

        # Check if you need to apply an aperture correction to the catalogue fluxes
        if 'aper_corr' in self.observed_catalogue.data.dtype.names:
            aper_corr = 10.**(-0.4*observation[0]['aper_corr'])
        else:
            aper_corr = 1.

        # Put observed photometry and its error in arrays
        obs_flux, obs_flux_err = self.observed_catalogue.extract_fluxes(self.filters, ID, key=self.key)

        # Sort wl_eff array
        wl_eff = np.array(self.filters.data['wl_eff'])
        sor = np.argsort(wl_eff)

        obs_flux, obs_flux_err, wl_eff = obs_flux[sor], obs_flux_err[sor], wl_eff[sor]

        obs_flux *= flux_factor
        obs_flux_err *= flux_factor

        # Consider only those bands with measurements!
        ok = np.where(obs_flux_err > 0.)[0]

        # Open the file containing BEAGLE results
        fits_file = os.path.join(BeagleDirectories.results_dir,
                str(ID)+'_BEAGLE.fits.gz')

        model_hdu = fits.open(fits_file)
        model_sed = model_hdu['marginal photometry']

        # Open the file containing the replicated data
        fits_file = os.path.join(BeagleDirectories.results_dir,
                BeagleDirectories.pypbeagle_data,
                str(ID)+'_BEAGLE_replic_data.fits.gz')

        replic_hdu = fits.open(fits_file)
        replic_data = replic_hdu[1]

        n_replicated = replic_data.data.field(0).size

        # the first column is the ID, so the number of bands is n-1
        n_bands = len(replic_data.columns.names)-1

        indices = replic_data.data['row_index']
        noiseless_flux = np.zeros((n_bands, indices.size))

        for i, band_name in enumerate((self.filters.data['label'][sor])):
            noiseless_flux[i, :] = model_sed.data['_'+band_name+'_'] / Jy * flux_factor

        # Compute the p-value band-by-band
        p_value_bands = np.zeros(n_bands)
        replic_fluxes = np.zeros((n_bands, n_replicated))
        for i, band_name in enumerate((self.filters.data['label'][sor])):
            
            if obs_flux_err[i] > 0.:
                obs_discr = (obs_flux[i].repeat(n_replicated)-noiseless_flux[i, :])**2 / obs_flux_err[i].repeat(n_replicated)**2
                repl_discr = (replic_data.data['_'+band_name+'_']*flux_factor - noiseless_flux[i, :])**2 / obs_flux_err[i].repeat(n_replicated)**2

                replic_fluxes[i,:] = replic_data.data['_'+band_name+'_'] * flux_factor
                p_value_bands[i] = 1. * np.count_nonzero((repl_discr >
                    obs_discr)) / n_replicated

        markers = np.array("o").repeat(n_bands)
        loc = np.where(p_value_bands <= p_value_lim)[0]
        markers[loc] = "o"
        print("p_value_bands: ", p_value_bands)

        ext_obs_flux = obs_flux.reshape(n_bands, 1).repeat(n_replicated, 1)
        ext_obs_flux_err = obs_flux_err.reshape(n_bands, 1).repeat(n_replicated, 1)

        obs_discr = np.sum((ext_obs_flux[ok,:]-noiseless_flux[ok,:])**2 / ext_obs_flux_err[ok,:]**2, axis=0)
        repl_discr = np.sum((replic_fluxes[ok,:]-noiseless_flux[ok,:])**2 / ext_obs_flux_err[ok,:]**2, axis=0)

        p_value = 1. * np.count_nonzero((repl_discr >
            obs_discr)) / n_replicated
        
        print("p_value: ", p_value)

        median_flux = np.zeros(n_bands)
        pdf_norm = np.zeros(n_bands)
        _max_y = np.zeros(n_bands)
        max_abs_flux = np.zeros(n_plot_x*n_plot_y)

        # Compute mean residual
        replic_fluxes = np.zeros((n_bands, n_replicated))
        mean_replic_fluxes = np.zeros(n_bands)
        for i, band_name in enumerate((self.filters.data['label'][sor])):
            mean_replic_fluxes[i] = np.mean(replic_data.data['_'+band_name+'_']) * flux_factor
            replic_fluxes[i,:] = replic_data.data['_'+band_name+'_'] * flux_factor

        mean_residual = (mean_replic_fluxes-obs_flux)/obs_flux_err 

        # Compute variance-covariance matrix of residual
        residual_fluxes = (replic_fluxes-obs_flux.reshape(n_bands,
            1).repeat(n_replicated, 1)) / obs_flux_err.reshape(n_bands,
                    1).repeat(n_replicated, 1)

        residual_covar = np.cov(residual_fluxes)
        print("residual_covar: ", residual_covar)

        # Plot the variance-covariance matrix of residuals
        if 'sns' in sys.modules:

            sns.set(style="white")
            labels = list()
            for lab in self.filters.data['label']:
                labels.append(lab.split('_')[-1])
                
            #d = pd.DataFrame(data=np.abs(residual_fluxes.T),
            #        columns=labels)

            # Compute the correlation matrix
            #corr = d.corr()

            # Generate a mask for the upper triangle
            mask = np.zeros_like(corr, dtype=np.bool)
            mask[np.triu_indices_from(mask)] = True

            # Set up the matplotlib figure
            fig, ax = plt.subplots()

            # Generate a custom diverging colormap
            cmap = sns.diverging_palette(220, 10, as_cmap=True)

            # Draw the heatmap with the mask and correct aspect ratio
            sns.heatmap(corr, mask=mask, cmap=cmap, vmax=0.4,
                    square=True, annot=True, fmt=".2f", annot_kws={"size": 10},
                    linewidths=.5, cbar_kws={"shrink": .85}, ax=ax)

            # Rotate by 45 deg the x and y ticks so they do not overlap
            plt.setp( ax.xaxis.get_majorticklabels(), rotation=45,
                    horizontalalignment='right' )

            plt.setp( ax.yaxis.get_majorticklabels(), rotation=45,
                    horizontalalignment='right' )

            name = prepare_plot_saving(str(ID)+'_BEAGLE_replic_data_phot_matrix.pdf')

            fig.savefig(name, dpi=None, facecolor='w', edgecolor='w',
                    orientation='portrait', papertype='a4', format="pdf",
                    transparent=False, bbox_inches="tight", pad_inches=0.1)

            fig.clear()
            plt.close(fig)

        # Select a random set of repliated data
        np.random.seed(seed=12345678)
        replic_data_rows = np.random.choice(n_replicated, size=n_plot_x*n_plot_y)    

        fig, axs = plt.subplots(n_plot_x, n_plot_y, sharex=True, sharey=True)
        fig.subplots_adjust(left=0.08, bottom=0.08, hspace=0, wspace=0)
        fontsize = 8
        axes_linewidth = 0.7

        ix = 0
        iy = 0
        for i, ax in enumerate(np.ravel(axs)):


#            kwargs = {'alpha':0.7}
#            (_, caps, _) = ax.errorbar(wl_eff, 
#                    obs_flux, 
#                    yerr=obs_flux_err, 
#                    ls=' ', 
#                    marker='o', 
#                    markersize=5, 
#                    color='orangered',
#                    markeredgewidth = 0.,
#                    elinewidth=1.0,
#                    capsize=2,
#                    **kwargs)
#
#            for cap in caps:
#                cap.set_color('orangered')
#                cap.set_markeredgewidth(1)
#
            temp_data = replic_data.data[replic_data_rows[i]]
            replic_fluxes = np.array(temp_data[1:]) * flux_factor

            diff_fluxes = (replic_fluxes-obs_flux) / obs_flux_err

            kwargs = {'alpha':0.4}
            unique_markers = np.unique(markers)
            if i != 0:
                for um in unique_markers:
                    mask = markers == um 

                    (_, caps, _) = ax.errorbar(wl_eff[mask], 
                            diff_fluxes[mask], 
                            ls=' ', 
                            marker=um, 
                            markersize=5, 
                            color='black',
                            markeredgewidth = 0.,
                            elinewidth=1.0,
                            capsize=2,
                            **kwargs)

                    for cap in caps:
                        cap.set_color('black')
                        cap.set_markeredgewidth(1)

            if i == 0:

                nXgrid = 1000
                kwargs = {'alpha':0.7}

                delta_wl = wl_eff[1:]-wl_eff[0:-1]
                delta_wl = np.concatenate(([delta_wl[0]], delta_wl))

                for j in range(n_bands):

                    residual = residual_fluxes[j,:]

                    # This function provides you with all the necessary info to draw violin plots
                    kde_pdf, pdf_norm, median_flux, x_plot, y_plot = prepare_violin_plot(residual)

                    w = 0.4 * delta_wl[j] / np.max(y_plot)

                    y_grid = np.full(len(x_plot), wl_eff[j])

                    _lim_y = kde_pdf(median_flux)/pdf_norm * w

                    ax.fill_betweenx(x_plot,
                            y_grid - y_plot*w,
                            y_grid + y_plot*w,
                            **kwargs
                            )

                    ax.plot( [wl_eff[j]-_lim_y, wl_eff[j]+_lim_y],
                            [median_flux, median_flux],
                            color = 'black',
                            linewidth = 0.2
                            )

            #min_flux[i] = np.min(np.array([replic_fluxes-obs_flux_err, obs_flux-obs_flux_err]))
            #max_flux[i] = np.max(np.array([replic_fluxes+obs_flux_err, obs_flux+obs_flux_err]))
            max_abs_flux[i] = np.max(np.abs(diff_fluxes))

        # Determine min and max values of y-axis
        yMax = np.max(max_abs_flux)
        yMin = -yMax
        dY = yMax-yMin
        yMax += dY * 0.1
        yMin -= dY * 0.1

        xMin = wl_eff[0]
        xMax = wl_eff[-1]
        dX = xMax-xMin
        xMax += dX * 0.1
        xMin -= dX * 0.1

        for ax in np.ravel(axs):        
            ax.set_ylim([yMin, yMax])
            ax.set_xlim([xMin, xMax])

            ax.tick_params(which='minor', axis='both',
                            length=2, width=axes_linewidth)

            ax.tick_params(which='major', axis='both',
                            length=3.5, width=axes_linewidth)

            x0, x1 = ax.get_xlim()
            if yMin < 0.: ax.plot( [x0,x1], [0.,0.], color='gray', lw=0.8 )


            for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                         ax.get_xticklabels() + ax.get_yticklabels()):
                item.set_fontsize(fontsize)


            for axis in ['top','bottom','left','right']:
              ax.spines[axis].set_linewidth(axes_linewidth)

            # Set better location of tick marks
            set_plot_ticks(ax, n_x=4, prune_x='both', prune_y='both')

        xlabel = "$\\lambda_\\textnormal{eff} / \\textnormal{\AA}$ (observed-frame)"
        #ylabel = "$f_{\\nu}/\\textnormal{nanoJy}$"
        ylabel = "$\\left(f_{\\nu}^\\textnormal{rep}-f_{\\nu}\\right) / \sigma$"

        fig.text(0.5, 0.02, xlabel, ha='center', fontsize=fontsize+1)
        fig.text(0.03, 0.5, ylabel, va='center', rotation='vertical', fontsize=fontsize+1)

        # Title of the plot is the object ID
        #plt.title(str(ID))

        # Location of printed text
        #x0, x1 = ax.get_xlim()
        #x = x0 + (x1-x0)*0.03
        #y0, y1 = ax.get_ylim()
        #y = y1 - (y1-y0)*0.10

        if print_text:

            # Print the evidence
            try:
                ax.text(x, y, "$\\log(Z)=" + "{:.2f}".format(self.logEvidence) + "$", fontsize=10 )
            except AttributeError:
                print("ciao")

            # Print the average reduced chi-square
            try:
                row = extract_row(self.PPC.data, ID, key=self.key)
                aver_chi_square = row['aver_chi_square']
                y = y1 - (y1-y0)*0.15
                ax.text(x, y, "$\\langle\chi^2\\rangle=" + "{:.2f}".format(aver_chi_square) + "$", fontsize=10 )
            except AttributeError:
                print("`PosteriorPredictiveChecks` not computed/loaded, hence " \
                "<chi^2> for the object `" + str(ID) + "` is not available")

            try:
                row = extract_row(self.PPC.data, ID, key=self.key)
                aver_red_chi_square = row['aver_red_chi_square']
                n_data = row['n_used_bands']
                y = y1 - (y1-y0)*0.20
                ax.text(x, y,
                        "$\\langle\chi^2/(\\textnormal{N}_\\textnormal{data}-1)\\rangle=" \
                        + "{:.2f}".format(aver_red_chi_square) + "\; \
                        (\\textnormal{N}_\\textnormal{data}=" + \
                        "{:d}".format(n_data) + ")" + "$", fontsize=10 )
            except AttributeError:
                print("`PosteriorPredictiveChecks` not computed/loaded, hence " \
                "<chi^2_red> for the object `" + str(ID) + "` is not available")


        #fig.tight_layout()

        name = prepare_plot_saving(plot_name)

        fig.savefig(name, dpi=None, facecolor='w', edgecolor='w',
                orientation='portrait', papertype='a4', format="pdf",
                transparent=False, bbox_inches="tight", pad_inches=0.1)

        plt.close(fig)

        model_hdu.close()
        replic_hdu.close()

##    def plot_residuals(self, residual_file_name=None, residual_plotname=None):
##
##        if not hasattr(self, 'observed_catalogue'):
##            except AttributeError:
##                    "An observed catalogue must be loaded before plotting the
##                    residual"
##
##        if not hasattr(self, 'beagle_summary_catalogue'):
##            except AttributeError:
##                    "A `beagle_summary_catalogue` must be loaded before plotting the
##                    residual"
##
##        self.residual = ResidualPhotometry()
##
##        try:
##            self.residual.load(self.residual_file_name)
##        except:
##            self.residual.compute(self.observed_catalogue,
##                self.beagle_summary_catalogue, self.self.filters.
##                cPickleName=self.residual_file_name)
                


