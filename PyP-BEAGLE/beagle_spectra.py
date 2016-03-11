import os
import logging
from scipy.integrate import simps, cumtrapz
from scipy.interpolate import interp1d
from bisect import bisect_left
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as plticker
import pandas as pd
# SEABORN creates by default plots with a filled background!!
#import seaborn as sns
from astropy.io import ascii
from astropy.io import fits

import sys
sys.path.append("../dependencies")
import WeightedKDE
#import FillBetweenStep

from beagle_utils import BeagleDirectories, prepare_plot_saving, set_plot_ticks, plot_exists
from beagle_filters import PhotometricFilters
from beagle_summary_catalogue import BeagleSummaryCatalogue
#from beagle_residual_photometry import ResidualPhotometry
from beagle_multinest_catalogue import MultiNestCatalogue
from beagle_posterior_predictive_checks import PosteriorPredictiveChecks


microJy = np.float32(1.E-23 * 1.E-06)
nanoJy = np.float32(1.E-23 * 1.E-09)

p_value_lim = 0.05

class ObservedSpectrum:

    def load(self, file_name):

        """ 
        Load an observed spectrum. It automatically
        detects, and loads, FITS or ASCII files depending on the suffix.

        Parameters
        ----------
        file_name : str
            Contains the file name of the spectrum.
        """

        if file_name.endswith(('fits', 'fit', 'FITS', 'FIT')):
            self.data = fits.open(file_name)[1].data
            self.columns = fits.open(file_name)[1].columns
        else:
            self.data = ascii.read(file_name, Reader=ascii.basic.CommentedHeader)


class Spectrum:

    def __init__(self):

        self.observed_spectrum = ObservedSpectrum()

        self.summary_catalogue = BeagleSummaryCatalogue()

        self.multinest_catalogue = MultiNestCatalogue()

        #self.residual = ResidualPhotometry()

        self.PPC = PosteriorPredictiveChecks()

    def plot_marginal(self, ID, max_interval=95.0,
            print_text=False, print_title=False, draw_steps=False, replot=False):    
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
        """

        # Name of the output plot
        plot_name = str(ID) + '_BEAGLE_marginal_SED_spec.pdf'

        # Check if the plot already exists
        if plot_exists(plot_name) and not replot:
            logging.warning('The plot "' + plot_name + '" already exists. \n Exiting the function.')
            return

        # The observed spectrum
        observation = self.observed_spectrum

        # Add to the error array the minimum relative error thet BEAGLE allows
        # one to add to the errors quoted in the catalogue
        #for i, err in enumerate(self.filters.data['flux_errcolName']):
        #    tmp_err = observation[0][err]
        #    if tmp_err > 0.:
        #        obs_flux_err[i] = observation[0][err]*aper_corr*self.filters.units / nanoJy
        #        obs_flux_err[i] = (np.sqrt( (obs_flux_err[i]/obs_flux[i])**2 +
        #                np.float32(self.filters.data['min_rel_err'][i])**2) *
        #                obs_flux[i])
        #    else:
        #        obs_flux_err[i] = tmp_err

        #ok = np.where(obs_flux_err > 0.)[0]

        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)

        # Open the file containing BEAGLE results
        fits_file = os.path.join(BeagleDirectories.results_dir,
                str(ID) + '_' + BeagleDirectories.suffix + '.fits.gz')

        hdulist = fits.open(fits_file)

        # Read the template wl array, and the 2D flux array
        model_wl = hdulist['marginal sed wl'].data['wl'][0,:]
        model_fluxes = hdulist['marginal sed'].data

        # Read the posterior probability
        probability = hdulist['posterior pdf'].data['probability']

        # Now compute for each wl bin the sorted fluxes, you will need this to
        # calculate the median and percentiles for each wl bin
        sort_indices = np.argsort(model_fluxes, axis=0)

        # Now it's time to compute the median (observed-frame) SED and its percentiles
        n_wl = model_fluxes.shape[1]

        median_flux = np.zeros(n_wl)
        lower_flux = np.zeros(n_wl)
        upper_flux = np.zeros(n_wl)

        for i in range(n_wl):

            # Compute the cumulative probability
            sort_ = sort_indices[:,i]
            cumul_pdf = cumtrapz(probability[sort_], model_fluxes[sort_,i], initial = 0.)
            cumul_pdf /= cumul_pdf[len(cumul_pdf)-1]
            
            # Get the interpolant of the cumulative probability
            f_interp = interp1d(cumul_pdf, model_fluxes[sort_,i])

            # The median corresponds to a cumulative probability = 0.5
            median_flux[i] = f_interp(0.5)

            # Compute the percentiles for the different credible regions
            lev = (1.-max_interval/100.)/2.
            lower_flux[i] = f_interp(lev)

            lev = 1.-(1.-max_interval/100.)/2.
            upper_flux[i] = f_interp(lev)
    
        # Set the plot limits from the minimum and maximum wl_eff
        dwl = model_wl[model_wl.size-1] - model_wl[0]
        wl_low = model_wl[0] - dwl*0.02
        wl_up = model_wl[model_wl.size-1] + dwl*0.02
            
        ax.set_xlim([1.,5.])
        ax.set_ylim([0., 1.2*np.amax(median_flux[np.isfinite(median_flux)])])

        # Define plotting styles
        ax.set_xlabel("$\lambda / \mu\\textnormal{m}$ (observed-frame)")
        ax.set_ylabel("$F_{\\lambda} / (\\textnormal{erg} \, \
                \\textnormal{s}^{-1} \, \\textnormal{cm}^{-2} \, \
                \\textnormal{\AA}^{-1})$")

        # Set better location of tick marks
        #set_plot_ticks(ax, n_x=5)

        kwargs = {'alpha':0.9}
        if ( draw_steps ):
            ax.step(observation.data['wl'],
                    observation.data['flux'],
                    where="mid",
                    color = "red",
                    linewidth = 1.50,
                    **kwargs
                    )
        else:
            ax.plot(observation.data['wl'],
                    observation.data['flux'],
                    color = "red",
                    linewidth = 2.50,
                    **kwargs
                    )

##        kwargs = { 'alpha': 0.4 }
##        if ( draw_steps ):
##            FillBetweenStep.fill_between_steps(ax,
##                    observation.data['wl'],
##                    observation.data['flux']-observation.data['flux_err'],
##                    observation.data['flux']+observation.data['flux_err'],
##                    step_where="mid",
##                    color = "red", 
##                    linewidth=0,
##                    interpolate=True,
##                    **kwargs)
##        else:
##            ax.fill_between(observation.data['wl'],
##                    observation.data['flux']-observation.data['flux_err'],
##                    observation.data['flux']+observation.data['flux_err'],
##                    facecolor = "red", 
##                    linewidth=0,
##                    interpolate=True,
##                    **kwargs)


        kwargs = { 'alpha': 0.7 }
        if ( draw_steps ):
            ax.step(model_wl/1.E+04,
                    median_flux,
                    where="mid",
                    color = "blue",
                    linewidth = 0.7,
                    **kwargs
                    )
        else:
            ax.plot(model_wl/1.E+04,
                    median_flux,
                    color = "blue",
                    linewidth = 2.00,
                    **kwargs
                    )

        kwargs = { 'alpha': 0.2 }
        if ( draw_steps ):
            FillBetweenStep.fill_between_steps(ax,
                    model_wl/1.E+04,
                    lower_flux[:], 
                    upper_flux[:],
                    step_where="mid",
                    color = "blue", 
                    linewidth=0,
                    **kwargs)
        else:
            ax.fill_between(model_wl/1.E+04,
                    lower_flux[:],
                    upper_flux[:],
                    facecolor = "blue", 
                    linewidth=0,
                    interpolate=True,
                    **kwargs)


        kwargs = { 'alpha': 0.8 }


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
                ax.text(x, y, "$\log(Z)=" + "{:.2f}".format(self.logEvidence) + "$", fontsize=10 )
            except AttributeError:
                print "ciao"

            # Print the average reduced chi-square
            try:
                aver_chi_square = self.PPC.data['aver_chi_square'][self.PPC.data['ID'] == ID]
                y = y1 - (y1-y0)*0.15
                ax.text(x, y, "$\langle\chi^2\\rangle=" + "{:.2f}".format(aver_chi_square) + "$", fontsize=10 )
            except AttributeError:
                print "`PosteriorPredictiveChecks` not computed/loaded, hence " \
                "<chi^2> for the object `" + str(ID) + "` is not available"

            try:
                aver_red_chi_square = self.PPC.data['aver_red_chi_square'][self.PPC.data['ID'] == ID]
                n_data = self.PPC.data['n_used_bands'][self.PPC.data['ID'] == ID]
                y = y1 - (y1-y0)*0.20
                ax.text(x, y,
                        "$\langle\chi^2/(\\textnormal{N}_\\textnormal{data}-1)\\rangle=" \
                        + "{:.2f}".format(aver_red_chi_square) + "\; \
                        (\\textnormal{N}_\\textnormal{data}=" + \
                        "{:d}".format(n_data) + ")" + "$", fontsize=10 )
            except AttributeError:
                print "`PosteriorPredictiveChecks` not computed/loaded, hence " \
                "<chi^2_red> for the object `" + str(ID) + "` is not available"

        if y0 < 0.: plt.plot( [x0,x1], [0.,0.], color='gray', lw=1.0 )

        name = prepare_plot_saving(plot_name)

        fig.savefig(name, dpi=None, facecolor='w', edgecolor='w',
                orientation='portrait', papertype='a4', format="pdf",
                transparent=False, bbox_inches="tight", pad_inches=0.1)
        plt.close(fig)

        hdulist.close()
