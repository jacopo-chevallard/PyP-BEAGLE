import os
import logging
import ConfigParser
from scipy.interpolate import interp1d
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
sys.path.append(os.path.join(os.environ['PYP_BEAGLE'], "dependencies"))
import WeightedKDE
#import FillBetweenStep

from beagle_utils import BeagleDirectories, prepare_plot_saving, set_plot_ticks, plot_exists
from beagle_filters import PhotometricFilters
from beagle_summary_catalogue import BeagleSummaryCatalogue
#from beagle_residual_photometry import ResidualPhotometry
from beagle_multinest_catalogue import MultiNestCatalogue
from beagle_posterior_predictive_checks import PosteriorPredictiveChecks
from beagle_mock_catalogue import BeagleMockCatalogue

TOKEN_SEP = ":"
microJy = np.float32(1.E-23 * 1.E-06)
nanoJy = np.float32(1.E-23 * 1.E-09)

p_value_lim = 0.05

class ObservedSpectrum(object):

    def __init__(self):

        self.description  = None

    def configure(self, param_file=None, config=None):

        if param_file is None:
            param_file = os.path.join(BeagleDirectories.results_dir, BeagleDirectories.param_file)

        if config is None:
            config = ConfigParser.SafeConfigParser()
            config.read(param_file)

        line = config.get('main', 'SPECTRUM FILE DESCRIPTION')

        self.description = { 
                "wl"        : {"colName" : None, "conversion" : None, "dispersion" : None, "type" : None}, 
                "flux"      : {"colName" : None, "conversion" : None}, 
                "fluxerr"   : {"colName" : None},  
                "sky"       : {"colName" : None},  
                "mask"      : {"colName" : None},  
                "redshift"  : {"keyword" : None},  
                "min_rel_err" : None,  
                }

        for key, value in self.description.iteritems():

            token = key + TOKEN_SEP

            if isinstance(value, dict):
                for in_key, in_value in value.iteritems():
                    in_token = token + in_key + TOKEN_SEP
                    if in_token in line:
                        self.description[key][in_key] = line.split(in_token)[1].split(" ")[0]

            else:
                if token in line:
                    self.description[key] = line.split(token)[1].split(" ")[0]
                

    def load(self, file_name):

        """ 
        Load an observed spectrum. It automatically
        detects, and loads, FITS or ASCII files depending on the suffix.

        Parameters
        ----------
        file_name : str
            Contains the file name of the spectrum.
        """

        if self.description is None:
            msg = ("The `ObservedSpectrum.description` must be set through "
                    "the `configure` method!")
            raise AttributeError(msg)


        hdu = fits.open(file_name)
        data = hdu[1].data

        self.data = dict()

        # Set the array containing the wavelength bins of the spectrum 
        try:
            self.data['wl'] = np.array(data[self.description["wl"]["colName"]])
        except:
            msg = ("`" + self.description["wl"]["colName"] + "` not found "
                    "in the current spectrum!")
            raise AttributeError(msg)

        if self.description["wl"]["conversion"] is not None:
            self.data['wl'] *= np.float(self.description["wl"]["conversion"])

        # Set the array containing the flux
        try:
            self.data['flux'] = np.array(data[self.description["flux"]["colName"]])
        except:
            msg = ("`" + self.description["flux"]["colName"] + "` not found "
                    "in the current spectrum!")
            raise AttributeError(msg)

        if self.description["flux"]["conversion"] is not None:
            self.data['flux'] *= np.float(self.description["flux"]["conversion"])

        # Set the array containing the flux error
        if self.description["fluxerr"]["colName"] is not None:
            try:
                self.data['fluxerr'] = np.array(data[self.description["fluxerr"]["colName"]])
            except:
                msg = ("`" + self.description["fluxerr"]["colName"] + "` not found "
                        "in the current spectrum!")
                raise AttributeError(msg)

        # Set the redshift
        if self.description["redshift"]["keyword"] is not None:
            self.data['redshift'] = np.float(hdu[1].header[self.description["redshift"]["keyword"]])

        hdu.close()

class Spectrum(object):

    def __init__(self, params_file):

        self.observed_spectrum = ObservedSpectrum()

        self.summary_catalogue = BeagleSummaryCatalogue()

        self.multinest_catalogue = MultiNestCatalogue()

        self.mock_catalogue = BeagleMockCatalogue(params_file)

        #self.residual = ResidualPhotometry()

        self.PPC = PosteriorPredictiveChecks()

    def plot_marginal(self, ID, 
            observation_name=None,
            max_interval=95.0,
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

        # If needed load the observed spectrum
        if observation_name is not None:
            print "\nobservation_name: ", observation_name, '\n'
            self.observed_spectrum.load(observation_name)
    
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
            # ******************************************************************
            # Here you must simply use `cumsum`, and not `cumtrapz` as in
            # beagle_utils.prepare_violin_plot, since the output of MultiNest are a set
            # of weights (which sum up to 1) associated to each set of parameters (the
            # `p_j` of equation 9 of Feroz+2009), and not a probability density (as the
            # MultiNest README would suggest).
            # ******************************************************************
            sort_ = sort_indices[:,i]
            cumul_pdf = np.cumsum(probability[sort_])
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
            ax.step(observation.data['wl']/1.E+04,
                    observation.data['flux'],
                    where="mid",
                    color = "red",
                    linewidth = 1.50,
                    **kwargs
                    )
        else:
            ax.plot(observation.data['wl']/1.E+04,
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
