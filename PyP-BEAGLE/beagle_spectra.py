from __future__ import absolute_import
from __future__ import print_function
import os
import logging
import six.moves.configparser
from collections import OrderedDict
from scipy.interpolate import interp1d
from bisect import bisect_left
import numpy as np
from itertools import tee

import json
import matplotlib.pyplot as plt
import matplotlib.ticker as plticker
import matplotlib.gridspec as gridspec
from matplotlib import rcParams
from matplotlib.axes import Axes
#import pandas as pd
# SEABORN creates by default plots with a filled background!!
#import seaborn as sns
from astropy.io import ascii
from astropy.io import fits

import sys
import pyp_beagle.dependencies.WeightedKDE 
from pyp_beagle.dependencies.walker_random_sampling import WalkerRandomSampling
import pyp_beagle.dependencies.autoscale
from pyp_beagle.dependencies import FillBetweenStep
import pyp_beagle.dependencies.set_shared_labels  as shLab

from .beagle_utils import BeagleDirectories, prepare_plot_saving, set_plot_ticks, plot_exists
from .beagle_filters import PhotometricFilters
from .beagle_summary_catalogue import BeagleSummaryCatalogue
#from beagle_residual_photometry import ResidualPhotometry
from .beagle_multinest_catalogue import MultiNestCatalogue
from .beagle_posterior_predictive_checks import PosteriorPredictiveChecks
from .beagle_mock_catalogue import BeagleMockCatalogue

# See here
# http://peak.telecommunity.com/DevCenter/PythonEggs#accessing-package-resources
# for an explanation on this approach to include data files
from pkg_resources import resource_stream
import six
from six.moves import zip
from six.moves import range

TOKEN_SEP = ":"
microJy = np.float32(1.E-23 * 1.E-06)
nanoJy = np.float32(1.E-23 * 1.E-09)

p_value_lim = 0.05


def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)

class ObservedSpectrum(object):

    def __init__(self):

        self.description  = None

    def configure(self, param_file=None, config=None):

        if param_file is None:
            param_file = os.path.join(BeagleDirectories.results_dir, BeagleDirectories.param_file)

        if config is None:
            config = six.moves.configparser.SafeConfigParser()
            config.read(param_file)

        line = config.get('main', 'SPECTRUM FILE DESCRIPTION')

        self.description = { 
                "wl"        : {"colName" : None, "conversion" : None, "dispersion" : None, "type" : None}, 
                "flux"      : {"colName" : None, "conversion" : None}, 
                "fluxerr"   : {"colName" : None},  
                "sky"       : {"colName" : None},  
                "mask"      : {"colName" : None},  
                "redshift"  : {"keyword" : "redshift"},  
                "min_rel_err" : None,  
                }

        for key, value in six.iteritems(self.description):

            token = key + TOKEN_SEP

            if isinstance(value, dict):
                for in_key, in_value in six.iteritems(value):
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
        self.data['redshift'] = None
        if self.description["redshift"]["keyword"] is not None:
            try:
                self.data['redshift'] = np.float(hdu[1].header[self.description["redshift"]["keyword"]])
            except:
                pass

        # Set the mask spectrum
        if self.description["mask"]["colName"] is not None:
            self.data['mask'] = np.array(data[self.description["mask"]["colName"]], dtype=bool)

        hdu.close()

class Spectrum(object):

    def __init__(self, params_file,
            **kwargs):

        self.inset_fontsize = BeagleDirectories.inset_fontsize_fraction * BeagleDirectories.fontsize

        self.observed_spectrum = ObservedSpectrum()

        self.multinest_catalogue = MultiNestCatalogue()

        self.mock_catalogue = kwargs.get('mock_catalogue')

        #self.residual = ResidualPhotometry()

        self.PPC = PosteriorPredictiveChecks()

        if kwargs.get('line_labels_json') is None:
            self.line_labels = json.load(resource_stream(__name__, 'files/emission_lines.json'), 
                    object_pairs_hook=OrderedDict)
        else:
            with open(kwargs.get('line_labels_json')) as f:    
                 self.line_labels = json.load(f, object_pairs_hook=OrderedDict)

        self.plot_line_labels = kwargs.get('plot_line_labels', False)

        self.resolution = kwargs.get('resolution')
    
        self.wl_range = kwargs.get('wl_range')

        self.wl_units = kwargs.get('wl_units', 'micron')

        self.wl_rest = kwargs.get('wl_rest', False)

        self.log_flux = kwargs.get('plot_log_flux', False)

        self.plot_full_SED = kwargs.get('plot_full_SED', False)

        self.show_residual = kwargs.get('show_residual', False)

        self.print_ID = kwargs.get('print_ID', False)

        self.draw_steps = kwargs.get('draw_steps', False)

        self.n_SED_to_plot = kwargs.get('n_SED_to_plot', 100)

        self.plot_suffix = kwargs.get('plot_suffix')

    def plot_marginal(self, ID, 
            observation_name=None,
            max_interval=95.0,
            print_text=False, 
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

        # Factor to convert angstrom to input units
        if self.wl_units == 'micron':
            wl_factor = 1.E+04
            xlabel = "\mu\\textnormal{m}"
        elif self.wl_units == 'nm':
            wl_factor = 1.E+01
            xlabel = "\\textnormal{nm}"
        elif self.wl_units == 'ang':
            wl_factor = 1.
            xlabel = "\\textnormal{\\AA}"
        else:
            raise ValueError("Wavelength units `" + self.wl_units + "` not recognised!")

        # If needed load the observed spectrum
        if observation_name is not None:
            self.observed_spectrum.load(observation_name)
    
        # Name of the output plot
        plot_suffix = '' 
        if self.plot_suffix is not None: plot_suffix = '_' + self.plot_suffix
        if self.wl_rest:
            plot_name = str(ID) + '_BEAGLE_marginal_SED_spec_rest_wl' + plot_suffix + '.pdf'
        else:
            plot_name = str(ID) + '_BEAGLE_marginal_SED_spec' + plot_suffix + '.pdf'

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
        residual_axs = list()
        if observation.data['redshift'] is not None:
            redshift = observation.data['redshift']
        else:
            _redshifts =  hdulist['galaxy properties'].data['redshift']
            _, _counts = np.unique(_redshifts, return_counts=True)
            if len(_counts) > 1:
                raise ValueError("The `redshift` of the object is not unique!")
            redshift = _redshifts[0]

        z1 = 1. + redshift

        if self.wl_rest:
            data_wl = observation.data['wl'] / z1
            data_flux = observation.data['flux'] * z1
            data_flux_err = observation.data['fluxerr'] * z1

            model_wl /= z1
            median_flux *= z1
            lower_flux *= z1
            upper_flux *= z1
        else:
            data_wl = observation.data['wl']
            data_flux = observation.data['flux']
            data_flux_err = observation.data['fluxerr']

        # Convert to correct wl units
        model_wl /= wl_factor
        data_wl /= wl_factor

        # Load the mask spectrum if exists
        mask_color = "grey"
        has_mask = False
        data_mask = np.ones(len(data_wl), dtype=bool)
        if "mask" in observation.data:
            data_mask = observation.data['mask']
            has_mask = True

        # Create masked versions of arrays
        slices = list()
        masked_array = np.ma.array(data_wl, mask=~data_mask)
        slices.append(np.ma.clump_unmasked(masked_array))
        slices.append(np.ma.clump_masked(masked_array))

        model_mask = np.ones(len(model_wl), dtype=bool)
        if "marginal sed mask" in hdulist:
            model_mask = np.array(hdulist['marginal sed mask'].data['mask'][0], dtype=bool)

        # Create masked versions of arrays
        slices_model = list()
        masked_array = np.ma.array(model_wl, mask=~model_mask)
        slices_model.append(np.ma.clump_unmasked(masked_array))
        slices_model.append(np.ma.clump_masked(masked_array))


        if self.show_residual:
            n_outer = 2
        else:
            n_outer = 1

        if self.wl_range is None:
            n_ranges = 1
        else:
            if len(self.wl_range) == 2:
                n_ranges = 1
            else:
                n_ranges = int(1.*len(self.wl_range)/2.)

        fig = plt.figure(figsize=(12,8))
        if self.show_residual:
            fig, axs_ = plt.subplots(n_outer, n_ranges, gridspec_kw = {'height_ratios':[3, 1]})
        else:
            fig, axs_ = plt.subplots(n_outer, n_ranges)

        fig.subplots_adjust(wspace=0.1, hspace=0.0)
        if self.show_residual:
            axs = axs_[0,:] 
        else:
            axs = axs_

        if self.show_residual:
            residual_axs = axs_[1,:]

        try:
            isinstance(axs[0], Axes)
        except:
            axs = [axs]
            residual_axs = [residual_axs]

        if n_ranges == 1:
            if self.wl_range is None:
                dwl = data_wl[-1]-data_wl[0]
                wl_low = data_wl[0] - dwl*0.025
                wl_up = data_wl[-1] + dwl*0.025
            else:
                wl_low, wl_up = self.wl_range[0], self.wl_range[1]

            axs[0].set_xlim([wl_low, wl_up])
            if self.show_residual:
                residual_axs[0].set_xlim([wl_low, wl_up])
        else:
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

            if self.show_residual:
                t = residual_axs[0].transAxes.transform([(0,0), (1,1)])
                t = residual_axs[0].get_figure().get_dpi() / (t[1,1] - t[0,1]) / 72
                d = 0.5*(rcParams['xtick.major.size']*t)

                wl_l = wl_ranges[0]
                for i, (ax_l, ax_r) in enumerate(pairwise(residual_axs)):
                    kwargs = dict(transform=ax_l.transAxes, color='k', clip_on=False)
                    ax_l.spines['right'].set_visible(False)
                    if not ax_l.spines['left'].get_visible():
                        ax_l.yaxis.set_ticks_position('none')
                    else:
                        ax_l.yaxis.tick_left()
                    ax_l.set_xlim(wl_l)
                    ax_l.plot((1-d, 1+d), (-d, +d), **kwargs)        
                    #ax_l.plot((1-d, 1+d), (1-d, 1+d), **kwargs)        

                    kwargs = dict(transform=ax_r.transAxes, color='k', clip_on=False)
                    ax_r.spines['left'].set_visible(False)
                    ax_r.yaxis.tick_right()
                    ax_r.tick_params(labelright='off') 
                    wl_r = wl_ranges[1+i]
                    ax_r.set_xlim(wl_r)
                    ax_r.plot((-d, +d), (-d, +d), **kwargs)        
                    #ax_r.plot((-d, +d), (1-d, 1+d), **kwargs)        
                    #ax_l.spines['top'].set_visible(False)
                    ax_l.spines['top'].set_color('none')
                    ax_l.xaxis.set_ticks_position('bottom')

                    ax_r.spines['top'].set_color('none')
                    ax_r.xaxis.set_ticks_position('bottom')
                    #ax_r.spines['top'].set_visible(False)

                    ax_r.patch.set_facecolor('None')
                    ax_l.patch.set_facecolor('None')

                    wl_l = wl_r

        which = 'both'
        for ax in axs:
            if self.log_flux:
                ax.set_yscale('log')
                which = 'x'

            # Set better location of tick marks
            if self.wl_range is not None:
                set_plot_ticks(ax, which=which, prune_x='both', n_x=3)
                for tick in ax.get_xticklabels():
                    tick.set_rotation(45)
                    if self.show_residual:
                        tick.set_ha('left')
                    else:
                        tick.set_ha('right')
                    tick.set_rotation_mode("anchor")
            else:
                set_plot_ticks(ax, which=which)

        if self.show_residual:
            for ax in axs:
                ax.tick_params(labelbottom='off')

            for ax in residual_axs:
                if self.wl_range is not None:
                    set_plot_ticks(ax, prune_x='both', n_x=3)
                    for tick in ax.get_xticklabels():
                        tick.set_rotation(45)

        # Define plotting styles
        xlabel = "$\\uplambda / " + xlabel + "$"
        if not self.wl_rest:
            xlabel = xlabel + " (observed frame)"

        if self.wl_range is not None:
            fig.text(0.5, -0.04, xlabel, ha='center')
        else:
            fig.text(0.5, 0.02, xlabel, ha='center')


        ylabel = "$F_{\\uplambda} / (\\textnormal{erg} \, \
                \\textnormal{s}^{-1} \, \\textnormal{cm}^{-2} \, \
                \\textnormal{\AA}^{-1})$"
        axs[0].set_ylabel(ylabel)

        if self.show_residual:
            ylabel = "$(F_{\\uplambda}-F_{\\uplambda}^\\textnormal{mod}) / F_{\\uplambda}$"
            residual_axs[0].set_ylabel(ylabel)

        # Title of the plot is the object ID
        if self.print_ID: 
            #fig.text(0.5, 0.95, str(ID).split('_')[0].strip(), ha='center', va='top')
            plt.suptitle(str(ID).split('_')[0].strip())

        alpha_line = 0.7
        alpha_fill = 0.3

        for ax in axs:

            colors = ["red", "salmon"]
            for slic, col in zip(slices, colors):
                for s in slic:
                    if (self.draw_steps):
                        ax.step(data_wl[s],
                                data_flux[s],
                                where="mid",
                                color=col,
                                linewidth=1.50,
                                alpha=alpha_line
                                )
                    else:
                        ax.plot(data_wl[s],
                                data_flux[s],
                                color=col,
                                linewidth=2.00,
                                alpha=alpha_line
                                )

                    if (self.draw_steps):
                        FillBetweenStep.fill_between_steps(ax,
                                data_wl[s],
                                data_flux[s]-data_flux_err[s],
                                data_flux[s]+data_flux_err[s],
                                step_where="mid",
                                color=col, 
                                linewidth=0,
                                interpolate=True,
                                alpha=alpha_fill
                                )
                    else:
                        ax.fill_between(data_wl[s],
                                data_flux[s]-data_flux_err[s],
                                data_flux[s]+data_flux_err[s],
                                facecolor=col, 
                                linewidth=0,
                                interpolate=True,
                                alpha=alpha_fill
                                )


            colors = ["blue", "deepskyblue"]
            for slic, col in zip(slices_model, colors):
                for s in slic:
                    if (self.draw_steps):
                        ax.step(model_wl[s],
                                median_flux[s],
                                where="mid",
                                color=col,
                                linewidth = 1.0,
                                alpha=alpha_line
                                )
                    else:
                        ax.plot(model_wl[s],
                                median_flux[s],
                                color=col,
                                linewidth = 1.5,
                                alpha=alpha_line
                                )

                    if (self.draw_steps):
                        FillBetweenStep.fill_between_steps(ax,
                                model_wl[s],
                                lower_flux[s], 
                                upper_flux[s],
                                step_where="mid",
                                color=col, 
                                linewidth=0,
                                alpha=alpha_fill
                                )
                    else:
                        ax.fill_between(model_wl[s],
                                lower_flux[s],
                                upper_flux[s],
                                facecolor=col, 
                                linewidth=0,
                                interpolate=True,
                                alpha=alpha_fill
                                )

            # Extract and plot full SED
            if 'full sed wl' in hdulist and self.plot_full_SED:
                indices = np.arange(len(probability))
                wrand = WalkerRandomSampling(probability, keys=indices)
                rand_indices = wrand.random(self.n_SED_to_plot)

                if self.wl_rest:
                    wl_obs = hdulist['full sed wl'].data['wl'][0,:]

                    for i in rand_indices:
                        SED = hdulist['full sed'].data[i,:]
                else:
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

        ylim = autoscale.get_autoscale_y(axs[0], ylog=self.log_flux)
        for ax in axs:
            yl = autoscale.get_autoscale_y(ax, ylog=self.log_flux, top_margin=0.2, bottom_margin=0.1)
            ylim = [min(ylim[0], yl[0]), max(ylim[1], yl[1])]

        for ax in axs:
            ax.set_ylim(ylim)

        if self.show_residual:
            data_flux_ = data_flux[data_wl==model_wl]
            data_mask_ = data_mask[data_wl==model_wl]
            data_flux_err_ = data_flux_err[data_wl==model_wl]
            residual = (data_flux_-median_flux)/data_flux_
            residual_err = (1./data_flux_ - (data_flux_-median_flux)/data_flux_**2) * data_flux_err_

            ymax = np.amax(abs(residual))
            ymax += 0.2*ymax
            ymax = 1.1

            colors = ["darkgreen"]
            indices = np.arange(len(data_wl))
            if has_mask:
                colors = ["darkgreen", mask_color]

            for ax in residual_axs:
                ax.set_ylim([-ymax, ymax])
                kwargs = {'alpha':0.7}

                ax.plot(ax.get_xlim(), [0.,0.],
                        color="black",
                        lw=0.7,
                        zorder=5,
                        **kwargs)

                for i, col in enumerate(colors):

                    if i==0:
                        ind = indices[data_mask_]
                    elif i==1:
                        ind = indices[~data_mask_]

                    ax.errorbar(model_wl[ind],
                            residual[ind],
                            yerr=residual_err[ind],
                            color = col,
                            ls=" ",
                            elinewidth = 0.5,
                            marker='o',
                            ms=3,
                            **kwargs
                            )

        for ax in axs:

            x0, x1 = ax.get_xlim()

            y0, y1 = ax.get_ylim()
            if self.log_flux:
                y0, y1 = np.log10(y0), np.log10(y1)
            dy = y1-y0

            # Label emission lines
            if self.plot_line_labels:
                for key, label in six.iteritems(self.line_labels):
                    if self.wl_rest:
                        x = label["wl"]/wl_factor
                    else:
                        x = label["wl"]/wl_factor * z1

                    if x < x0 or x > x1:
                        continue

                    if self.resolution is not None:
                        if abs(x-prev_x) < (x/self.resolution):
                            continue

                    direction = "up"
                    if 'direction' in label:
                        direction = label['direction']

                    color = "black"
                    if 'color' in label:
                        color = label['color']

                    i1 = bisect_left(data_wl, x-2)
                    i2 = bisect_left(data_wl, x+2.)
                    if direction == "up":
                        mm = np.amax(data_flux[i1:i2+1])
                        y = mm+dy*0.1
                        ydots = mm+dy*0.07
                    elif direction == "down":
                        mm = np.amin(data_flux[i1:i2+1])
                        y = mm-dy*0.1
                        ydots = mm-dy*0.07

                    ha = 'center'
                    if 'ha' in label:
                        ha = label['ha']

                    ax.text(x, y, 
                            label["label"], 
                            fontsize=self.inset_fontsize, 
                            rotation=45,
                            color=color,
                            va='center',
                            ha=ha)

                    y0 = mm
                    ax.plot([x,x], [y0,ydots], 
                            ls=":",
                            color=color,
                            lw=1.0,
                            zorder=0)

            if print_text:

                # Print the evidence
                try:
                    ax.text(x, y, "$\\log(Z)=" + "{:.2f}".format(self.logEvidence) + "$", fontsize=self.inset_fontsize)
                except AttributeError:
                    print("ciao")

                # Print the average reduced chi-square
                try:
                    aver_chi_square = self.PPC.data['aver_chi_square'][self.PPC.data['ID'] == ID]
                    y = y1 - (y1-y0)*0.15
                    ax.text(x, y, "$\\langle\chi^2\\rangle=" + "{:.2f}".format(aver_chi_square) + "$", fontsize=self.inset_fontsize)
                except AttributeError:
                    print("`PosteriorPredictiveChecks` not computed/loaded, hence " \
                    "<chi^2> for the object `" + str(ID) + "` is not available")

                try:
                    aver_red_chi_square = self.PPC.data['aver_red_chi_square'][self.PPC.data['ID'] == ID]
                    n_data = self.PPC.data['n_used_bands'][self.PPC.data['ID'] == ID]
                    y = y1 - (y1-y0)*0.20
                    ax.text(x, y,
                            "$\\langle\chi^2/(\\textnormal{N}_\\textnormal{data}-1)\\rangle=" \
                            + "{:.2f}".format(aver_red_chi_square) + "\; \
                            (\\textnormal{N}_\\textnormal{data}=" + \
                            "{:d}".format(n_data) + ")" + "$", fontsize=self.inset_fontsize)
                except AttributeError:
                    print("`PosteriorPredictiveChecks` not computed/loaded, hence " \
                    "<chi^2_red> for the object `" + str(ID) + "` is not available")

            if y0 < 0.: plt.plot( [x0,x1], [0.,0.], color='gray', lw=1.0 )

        #plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
        #plt.xlabel(
        #plt.ylabel("$F_{\\lambda} / (\\textnormal{erg} \, \
        #        \\textnormal{s}^{-1} \, \\textnormal{cm}^{-2} \, \
        #        \\textnormal{\AA}^{-1})$", va='center', rotation='vertical')

        #fig.text(0.5, 0.02, , ha='center')
        #fig.text(0.04, 0.5, 
        #        \\textnormal{s}^{-1} \, \\textnormal{cm}^{-2} \, \
        #        \\textnormal{\AA}^{-1})$", va='center', rotation='vertical')


        #plt.tight_layout()

        name = prepare_plot_saving(plot_name)

        fig.savefig(name, dpi=None, facecolor='w', edgecolor='w',
                orientation='portrait', papertype='a4', format="pdf",
                transparent=False, bbox_inches="tight", pad_inches=0.1)
        plt.close(fig)

        hdulist.close()
