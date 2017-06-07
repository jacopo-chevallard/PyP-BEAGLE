import os
import logging
import ConfigParser
from collections import OrderedDict
from scipy.interpolate import interp1d
from bisect import bisect_left
import numpy as np
from itertools import tee, izip

import json
import matplotlib.pyplot as plt
import matplotlib.ticker as plticker
from matplotlib import rcParams
#import pandas as pd
# SEABORN creates by default plots with a filled background!!
#import seaborn as sns
from astropy.io import ascii
from astropy.io import fits

import sys
import dependencies.WeightedKDE as WeightedKDE
import dependencies.autoscale as autoscale
from dependencies.walker_random_sampling import WalkerRandomSampling
#import FillBetweenStep

from beagle_utils import BeagleDirectories, prepare_plot_saving, set_plot_ticks, plot_exists
from beagle_filters import PhotometricFilters
from beagle_summary_catalogue import BeagleSummaryCatalogue
#from beagle_residual_photometry import ResidualPhotometry
from beagle_multinest_catalogue import MultiNestCatalogue
from beagle_posterior_predictive_checks import PosteriorPredictiveChecks
from beagle_mock_catalogue import BeagleMockCatalogue

# See here
# http://peak.telecommunity.com/DevCenter/PythonEggs#accessing-package-resources
# for an explanation on this approach to include data files
from pkg_resources import resource_stream

TOKEN_SEP = ":"
microJy = np.float32(1.E-23 * 1.E-06)
nanoJy = np.float32(1.E-23 * 1.E-09)

p_value_lim = 0.05


def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return izip(a, b)

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


        hdu = fits.open(os.path.expandvars(file_name))
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

    def __init__(self, params_file,
            line_labels_json=None,
            plot_line_labels=False,
            resolution=None,
            mock_catalogue=None,
            wl_range=None,
            plot_full_SED=False,
            n_SED_to_plot=100):

        self.observed_spectrum = ObservedSpectrum()

        self.multinest_catalogue = MultiNestCatalogue()

        self.mock_catalogue = mock_catalogue

        #self.residual = ResidualPhotometry()

        self.PPC = PosteriorPredictiveChecks()

        if line_labels_json == None:
            self.line_labels = json.load(resource_stream(__name__, 'files/emission_lines.json'), 
                    object_pairs_hook=OrderedDict)
        else:
            with open(line_labels_json) as f:    
                 self.line_labels = json.load(f, object_pairs_hook=OrderedDict)

        self.plot_line_labels = plot_line_labels

        self.resolution = resolution
    
        self.wl_range = wl_range

        self.plot_full_SED = plot_full_SED

        self.n_SED_to_plot = n_SED_to_plot

    def plot_marginal(self, ID, 
            observation_name=None,
            max_interval=95.0,
            print_text=False, 
            print_title=False, 
            draw_steps=False, 
            replot=False):    
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
        axs = list()
        if self.wl_range is None:
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1)
            dwl = observation.data['wl'][-1]-observation.data['wl'][0]
            wl_low = observation.data['wl'][0] - dwl*0.025
            wl_up = observation.data['wl'][-1] + dwl*0.025
            ax.set_xlim([wl_low/1.E+04, wl_up/1.E+04])
            axs.append(ax)
        else:
            if len(self.wl_range) == 2:
                fig = plt.figure()
                ax = fig.add_subplot(1, 1, 1)
                ax.set_xlim(self.wl_range)
                axs.append(ax)
            else:
                n_ranges = int(1.*len(self.wl_range)/2.)
                fig, axs = plt.subplots(1, n_ranges, sharey=True)
                fig.subplots_adjust(wspace=.1)

                # how big to make the diagonal lines in axes coordinates
                # converting "points" to axes coordinates: 
                # https://stackoverflow.com/a/33638091
                t = axs[0].transAxes.transform([(0,0), (1,1)])
                t = axs[0].get_figure().get_dpi() / (t[1,1] - t[0,1]) / 72
                d = 0.5*(rcParams['xtick.major.size']*t)

                wl_ranges = [(self.wl_range[2*i], self.wl_range[2*i+1]) for i in range(n_ranges)]

                wl_l = wl_ranges[0]
                for i, (ax_l, ax_r) in enumerate(pairwise(axs)):
                    kwargs = dict(transform=ax_l.transAxes, color='k', clip_on=False)
                    ax_l.spines['right'].set_visible(False)
                    if not ax_l.spines['left'].get_visible():
                        ax_l.yaxis.set_ticks_position('none')
                    else:
                        ax_l.yaxis.tick_left()
                    ax_l.set_xlim(wl_l)
                    ax_l.plot((1-d, 1+d), (-d, +d), **kwargs)        
                    ax_l.plot((1-d, 1+d), (1-d, 1+d), **kwargs)        

                    kwargs = dict(transform=ax_r.transAxes, color='k', clip_on=False)
                    ax_r.spines['left'].set_visible(False)
                    ax_r.yaxis.tick_right()
                    ax_r.tick_params(labelright='off') 
                    wl_r = wl_ranges[1+i]
                    ax_r.set_xlim(wl_r)
                    ax_r.plot((-d, +d), (-d, +d), **kwargs)        
                    ax_r.plot((-d, +d), (1-d, 1+d), **kwargs)        

                    wl_l = wl_r

        for ax in axs:
            ax.set_ylim([0., 1.2*np.amax(median_flux[np.isfinite(median_flux)])])

            # Set better location of tick marks
            if len(axs) > 1:
                set_plot_ticks(ax, prune_x='both')
            else:
                set_plot_ticks(ax)

        # Define plotting styles
        fig.text(0.5, 0.02, "$\lambda / \mu\\textnormal{m}$ (observed-frame)", ha='center')
        fig.text(0.04, 0.5, "$F_{\\lambda} / (\\textnormal{erg} \, \
                \\textnormal{s}^{-1} \, \\textnormal{cm}^{-2} \, \
                \\textnormal{\AA}^{-1})$", va='center', rotation='vertical')


        for ax in axs:

            kwargs = {'alpha':0.7}
            if (draw_steps):
                ax.step(observation.data['wl']/1.E+04,
                        observation.data['flux'],
                        where="mid",
                        color = "red",
                        linewidth = 2.50,
                        **kwargs
                        )
            else:
                ax.plot(observation.data['wl']/1.E+04,
                        observation.data['flux'],
                        color = "red",
                        linewidth = 2.50,
                        **kwargs
                        )

            kwargs = { 'alpha': 0.3 }
            if (draw_steps):
                FillBetweenStep.fill_between_steps(ax,
                        observation.data['wl']/1.E+04,
                        observation.data['flux']-observation.data['fluxerr'],
                        observation.data['flux']+observation.data['fluxerr'],
                        step_where="mid",
                        color = "red", 
                        linewidth=0,
                        interpolate=True,
                        **kwargs)
            else:
                ax.fill_between(observation.data['wl']/1.E+04,
                        observation.data['flux']-observation.data['fluxerr'],
                        observation.data['flux']+observation.data['fluxerr'],
                        facecolor = "red", 
                        linewidth=0,
                        interpolate=True,
                        **kwargs)


            kwargs = { 'alpha': 0.7 }
            if (draw_steps):
                ax.step(model_wl/1.E+04,
                        median_flux,
                        where="mid",
                        color = "blue",
                        linewidth = 1.5,
                        **kwargs
                        )
            else:
                ax.plot(model_wl/1.E+04,
                        median_flux,
                        color = "blue",
                        linewidth = 1.5,
                        **kwargs
                        )

            kwargs = { 'alpha': 0.3 }
            if (draw_steps):
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

            # Extract and plot full SED
            if 'full sed wl' in hdulist and self.plot_full_SED:
                indices = np.arange(len(probability))
                wrand = WalkerRandomSampling(probability, keys=indices)
                rand_indices = wrand.random(self.n_SED_to_plot)

                z1 = (1.+self.observed_spectrum.data['redshift'])
                wl_obs = hdulist['full sed wl'].data['wl'][0,:] * z1

                for i in rand_indices:
                    SED = hdulist['full sed'].data[i,:] / z1

                    ax.plot(wl_obs, 
                            flux_obs,
                            color="black",
                            ls="-",
                            lw=0.5,
                            alpha=0.5)


            kwargs = { 'alpha': 0.8 }

            #autoscale.autoscale_y(ax)

        # Title of the plot is the object ID
        if print_title: plt.title(str(ID))

        for ax in axs:

            # Location of printed text
            x0, x1 = ax.get_xlim()
            x = x0 + (x1-x0)*0.03
            y0, y1 = ax.get_ylim()
            ax.set_ylim([0.,y1])
            y = y1 - (y1-y0)*0.10


            # Label emission lines
            if self.plot_line_labels:
                x0, x1 = ax.get_xlim()
                n = 3
                i = 0
                prev_x = 0.
                for key, label in self.line_labels.iteritems():
                    x = label["wl"]/1.E+04 * (1.+self.observed_spectrum.data['redshift'])
                    if x < x0 or x > x1:
                        continue
                    if self.resolution is not None:
                        if abs(x-prev_x) < (x/self.resolution):
                            continue
                    if i >= n:
                        i=0
                    prev_x = x
                    y0, y1 = ax.get_ylim() ; dy = y1-y0
                    y = y1 - dy*0.05 - dy*0.025*i
                    ax.text(x, y, 
                            label["label"], 
                            fontsize=7, 
                            rotation=45,
                            ha='center')

                    y -= dy*0.025
                    ax.plot([x,x], [y,0.], 
                            ls="--",
                            lw=0.2,
                            color="black",
                            zorder=0)
                    i+=1


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
