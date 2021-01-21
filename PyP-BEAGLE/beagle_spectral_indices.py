from __future__ import absolute_import
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

from .significant_digits import to_precision

from .beagle_utils import BeagleDirectories, prepare_plot_saving, set_plot_ticks, plot_exists, \
        prepare_violin_plot, extract_row

from .beagle_observed_catalogue import ObservedCatalogue

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


_TOKENS = ['lum', 'lumErr', 'ew', 'ewErr']
_COL_NAME = 'colName'
_LABEL = 'label'
_SEP = ':'

class SpectralIndicesCatalogue(ObservedCatalogue):

    def configure(self, file_name):

        self.line_config = OrderedDict()

        for line in open(file_name):
            if line.strip() and not line.startswith("#"):
                _line = line

                _key = _LABEL + _SEP
                if _key in line:
                    _label = line.split(_key)[1].split()[0]
                    self.line_config[_label] = OrderedDict()

                    for t in _TOKENS:
                        _key = t + _SEP + _COL_NAME + _SEP
                        if _key in line:
                            _value = line.split(_key)[1].split()[0]
                            self.line_config[_label][t] = _value



class SpectralIndices(object):

    def __init__(self, **kwargs):

        with open(kwargs.get('line_labels_json')) as f:    
            self.line_list = json.load(f, object_pairs_hook=OrderedDict)

        self.inset_fontsize = BeagleDirectories.inset_fontsize_fraction * BeagleDirectories.fontsize

        self.observed_catalogue = SpectralIndicesCatalogue()

        self.key = kwargs.get('ID_key') 

        self.plot_log_flux = kwargs.get('plot_log_flux')

        self.print_values = kwargs.get('print_line_values')

        self.max_interval = kwargs.get('max_interval')


    def plot_line_fluxes(self, 
            ID, replot=False, 
            title=False, letter=None, signif_digits=1):

        suffix = ""
        #if self.print_values:
        #    suffix = "_printed_values"

        plot_name = str(ID) + '_BEAGLE_spectral_indices' + suffix + ".pdf"

        # Check if the plot already exists
        if plot_exists(plot_name) and not replot:
            logging.warning('The plot "' + plot_name + '" already exists. \n Exiting the function.')
            return

        # From the (previously loaded) observed catalogue select the row
        # corresponding to the input ID
        observation = extract_row(self.observed_catalogue.data, ID, key=self.key)

        fig = plt.figure(figsize=(12, 3))
        ax = fig.add_subplot(1, 1, 1)

        # Open the file containing BEAGLE results
        fits_file = os.path.join(BeagleDirectories.results_dir,
                str(ID) + '_' + BeagleDirectories.suffix + '.fits.gz')

        hdulist = fits.open(fits_file)

        # Read model fluxes
        model_fluxes = hdulist['spectral indices'].data
        n_lines = len(self.line_list)

        # Read the posterior probability
        probability = hdulist['posterior pdf'].data['probability']

        width = 0.5

        minY_values = np.zeros(n_lines) ; maxY_values = np.zeros(n_lines)
        _observed_fluxes = np.zeros(n_lines) ; _observed_flux_errors = np.zeros(n_lines)
        _model_fluxes = np.zeros(n_lines)
        for i, (key, value) in enumerate(six.iteritems(self.line_list)):

            X = i+1
            _line_conf = self.observed_catalogue.line_config[key]
            if 'lum' in _line_conf:
                _col_name = _line_conf['lum']
                _err_col_name = _line_conf['lumErr']
            elif 'ew' in _line_conf:
                _col_name = _line_conf['ew']
                _err_col_name = _line_conf['ewErr']

            _observed_flux = observation[_col_name]
            _observed_fluxes[i] = _observed_flux 

            _observed_flux_err = observation[_err_col_name]
            _observed_flux_errors[i] = _observed_flux_err

            _label = value["label"]

            # This function provides you with all the necessary info to draw violin plots
            _model_flux = model_fluxes[key]
            kde_pdf, pdf_norm, median_flux, x_plot, y_plot = prepare_violin_plot(_model_flux, weights=probability) 
            _model_fluxes[i] = median_flux

            _max_y = np.max(y_plot)
            w = width / _max_y

            y_grid = np.full(len(x_plot), X)
            _lim_y = kde_pdf(median_flux)/pdf_norm * w

            kwargs = {'alpha':0.8}
            ax.errorbar(X,
                    _observed_flux, 
                    yerr=_observed_flux_err,
                    color="dodgerblue",
                    ls=" ",
                    marker=" ",
                    zorder=5,
                    **kwargs)

            ax.plot(X,
                    _observed_flux, 
                    color="dodgerblue",
                    ls=" ",
                    marker="D",
                    markeredgewidth=0.,
                    markersize = 8,
                    zorder=3,
                    **kwargs)

            kwargs = {'color':'tomato', 'alpha':0.7, 'edgecolor':'black', 'linewidth':0.2}
            ax.fill_betweenx(x_plot,
                    y_grid - y_plot*w,
                    y_grid + y_plot*w,
                    zorder=2,
                    **kwargs
                    )

            ax.plot([X-_lim_y, X+_lim_y],
                    [median_flux, median_flux],
                    color = 'black',
                    zorder=2,
                    linewidth = 0.2
                    )

            ax.plot(X,
                    median_flux,
                    color = 'black',
                    marker = "o",
                    markersize = 5,
                    zorder=4,
                    alpha = 0.7
                    )

            _all = np.concatenate((_observed_flux-_observed_flux_err, x_plot))
            _min = np.amin(_all[_all > 0.])
            minY_values[i] = _min

            _all = np.concatenate((_observed_flux+_observed_flux_err, x_plot))
            _max = np.amax(_all[_all > 0.])
            maxY_values[i] = _max

        minY = np.amin(minY_values)
        maxY = np.max(maxY_values)
        if self.print_values:
            _factor = 0.15
        else:
            _factor = 0.10

        if self.plot_log_flux:
            dY = np.log10(maxY)-np.log10(minY)
            minY = 10.**(np.log10(minY)-_factor*dY)
            maxY = 10.**(np.log10(maxY)+_factor*dY)
        else:
            dY = maxY-minY
            minY -= _factor*dY
            maxY += _factor*dY
            
        ax.set_ylim(minY, maxY) 
        
        # Set better location of tick marks
        set_plot_ticks(ax, n_x=5)


        if self.plot_log_flux:
            ax.set_ylabel("$\\log(\\textnormal{F}/\\textnormal{erg} \; \
                \\textnormal{s}^{-1} \, \\textnormal{cm}^{-2})$")
        else:
            ax.set_ylabel("$\\textnormal{F}/\\textnormal{erg} \; \
                \\textnormal{s}^{-1} \, \\textnormal{cm}^{-2}$")


        plt.tick_params(
            axis='x',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            top='off',      # ticks along the bottom edge are off
            bottom='off')         # ticks along the top edge are off

        xticks = list(range(1, len(self.line_list)+1))
        ax.set_xticks(xticks)

        ticklabels = list()
        alpha = 0.6
        minY, maxY = ax.get_ylim() ; dY = maxY - minY
        _factor = 0.03  ; _init_fact = 0.04
        j = 1
        _fontsize = min(self.inset_fontsize, self.inset_fontsize/(0.075*n_lines))
        for i, (line_key, line_value) in enumerate(six.iteritems(self.line_list)):
            X = i+1
            ticklabels.append(line_value["label"])
            _variab_fact = _init_fact * j
            if self.plot_log_flux:
                dY = np.log10(maxY) - np.log10(minY)
                Y = 10.**(np.log10(minY_values[i]) - _factor*dY)
                Y_t0 = 10.**(np.log10(maxY) - _variab_fact*dY)
                Y_t1 = 10.**(np.log10(maxY) - _variab_fact*dY - 1.4*_init_fact*dY)
            else:
                Y = minY_values[i] - _factor*dY
                Y_t0 = maxY - _variab_fact*dY
                Y_t1 = maxY - _variab_fact*dY - 1.4*_init_fact*dY

            ax.plot([X,X], [minY, Y],
                    color="black",
                    ls=":",
                    zorder=2,
                    alpha=alpha)

            if self.print_values:
                _val = _observed_fluxes[i]
                if not _val > 0.:
                    continue
                _val_err = _observed_flux_errors[i]
                _n = int(np.floor(np.log10(_val)))
                _norm = 10.**_n
                #_text = '$' + to_precision(_val/_norm, signif_digits+1) + '\\pm' + to_precision(_val_err/_norm, signif_digits) + '\\; 10^{' + str(_n) + '}$'
                _text = '$' + to_precision(_val/_norm, signif_digits+2) + '\\pm' + to_precision(_val_err/_norm, signif_digits) + '$'
                ax.text(X, Y_t0,
                        _text,
                        horizontalalignment='center',
                        verticalalignment='center',
                        color="black", 
                        fontsize=_fontsize)

                _val = _model_fluxes[i]
                _text = '$' + to_precision(_val/_norm, signif_digits+2) + '$'
                ax.text(X, Y_t1,
                        _text,
                        horizontalalignment='center',
                        verticalalignment='center',
                        color="black", 
                        fontsize=_fontsize)

            j += 1
            if j%4 == 0:
                j = 1

        ax.set_xticklabels(ticklabels, rotation=45, rotation_mode="anchor", ha='right')
        #zed = [tick.label.set_fontsize(fontsize*factor) for tick in ax.yaxis.get_major_ticks()]

        if self.plot_log_flux:
            ax.set_yscale('log')

        if title:
            plt.title(title)

        if letter is not None:
            ax.text(0.010, 0.940, '('+letter+')',
                    horizontalalignment='left',
                    verticalalignment='center',
                    color="black", 
                    transform=ax.transAxes)

        name = prepare_plot_saving(plot_name)

        fig.savefig(name, dpi=None, facecolor='w', edgecolor='w',
                orientation='portrait', papertype='a4', format="pdf",
                transparent=False, bbox_inches="tight", pad_inches=0.1)

        plt.close(fig)
        hdulist.close()


