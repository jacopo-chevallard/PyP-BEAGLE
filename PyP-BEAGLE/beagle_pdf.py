import logging
import numpy as np
import os
import json
from matplotlib import rc
from matplotlib.colors import colorConverter
from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt
from collections import OrderedDict
from astropy.io import fits

from getdist import plots, MCSamples

from beagle_utils import BeagleDirectories, prepare_plot_saving, plot_exists

class PDF(object):

    def __init__(self, params_file):

        # Names of parameters, used to label the axes, and whether they are log
        # or not
        with open(params_file) as f:    
            # The use of the "OrderedDict" ensures that the order of the
            # entries in the dictionary reflects the order in the file
            self.adjust_params = json.load(f, object_pairs_hook=OrderedDict)

    def plot_triangle(self, ID, 
            mock_file_name=None,
            params_to_plot=None, 
            suffix=None, 
            replot=False, 
            M_star=False, 
            show=False):
        """ 
        Draw a "triangle plot" with the 1D and 2D posterior probability

        Parameters
        ----------
        ID : str or int
            Object ID


        M_star : bool, optional
            If set, the routine will plot the mass currenly locked into stars
            instead of the mass of star formed (i.e., the plotted mass will
            accout for the return fraction)

        """ 
        # NB: you changed the getdist/plot.py _set_locator function
        # replacing line 1172-1176
        #if x and (abs(xmax - xmin) < 0.01 or max(abs(xmin), abs(xmax)) >= 1000):
        #    axis.set_major_locator(plt.MaxNLocator(self.settings.subplot_size_inch / 2 + 3, prune=prune))
        #else:
        #    axis.set_major_locator(plt.MaxNLocator(self.settings.subplot_size_inch / 2 + 4, prune=prune))
        # with
        #if x and (abs(xmax - xmin) < 0.01 or max(abs(xmin), abs(xmax)) >= 1000):
        #    axis.set_major_locator(plt.MaxNLocator(self.settings.subplot_size_inch / 2 + 2, prune=prune))
        #else:
        #    axis.set_major_locator(plt.MaxNLocator(self.settings.subplot_size_inch / 2 + 3, prune=prune))
        # to have a fewer number of tick marks, and hence less crowded triangle plots

        # for "cosmetics" reasons you also changed the getdist/plot.py setWithSubplotSize function
        # replacing line 166-167
        # self.lab_fontsize = 7 + 2 * self.subplot_size_inch
        # self.axes_fontsize = 4 + 2 * self.subplot_size_inch
        # with
        # self.lab_fontsize = 7 + 4 * self.subplot_size_inch
        # self.axes_fontsize = 4 + 4 * self.subplot_size_inch
        # to have larger, hence more readable, label and axes font sizes
    
        # Name of the output plot
        if suffix is None:
            plot_name = str(ID)+'_BEAGLE_triangle.pdf'
        else:
            plot_name = str(ID)+'_BEAGLE_triangle_' + suffix + '.pdf'

        # Check if the plot already exists
        if plot_exists(plot_name) and not replot and not show:
            logging.warning('The plot "' + plot_name + '" already exists. \n Exiting the function.')
            return

        # Check if a file name containing the "true" value of the adjustable parameters has been provided
        mock_param_values = None
        if mock_file_name is not None:
            mock_param_values = list()
            hdulist = fits.open(mock_file_name)

            # This cycles over all keys containing the parameter names
            for key, value in self.adjust_params.iteritems():
                if "extName" in value:
                    extName = value["extName"]
                else:
                    extName = "POSTERIOR PDF"

                d = { "extName" : extName, "colName" : key}
                log = {"log" : False}
                log.update(value)
                # This cycles over the dictioary items for one parameter
                for in_key, in_value in value.iteritems():
                    if in_key == "mock":
                        # This will merge the default dictionary `d` with the
                        # one found in the json file
                        d.update(in_value)

                val = hdulist[d["extName"]].data[d["colName"]]
                if log["log"]:
                    val = np.log10(val)
                mock_param_values.append(val)
            hdulist.close()


        fits_file = os.path.join(BeagleDirectories.results_dir,
                str(ID) + '_' + BeagleDirectories.suffix + '.fits.gz')

        hdulist = fits.open(fits_file)

        param_values = OrderedDict()
        for key, value in self.adjust_params.iteritems():
            extName = "POSTERIOR PDF"
            if "extName" in value:
                extName = value["extName"]

            param_values[key] = hdulist[extName].data[key]

        probability = hdulist['posterior pdf'].data['probability']

        n_rows = probability.size

       # ParamsToPlot = ['mass', 'redshift', 'tauV_eff', 'metallicity', 'specific_sfr', 'tau']

        # By default you plot all parameters
        if params_to_plot is None:
            _params_to_plot = list()
            for key, value in self.adjust_params.iteritems():
                _params_to_plot.append(key)
        else: 
            _params_to_plot = params_to_plot

        # Here you check whether you want to plot the mass currently locked
        # into stars or not (i.e. accounting for the return fraction as well)
        if M_star and 'mass' in _params_to_plot:
            param_values['mass'][:] = np.log10(hdulist['galaxy properties'].data['M_star'][:])

        nParamsToPlot = len(_params_to_plot)

        names = list()
        labels = list()
        ranges = dict()
        samps = np.zeros((n_rows, nParamsToPlot))
        keys = list()

        j = 0
        for key, par in self.adjust_params.iteritems():
            keys.append(key)
            for par_name in _params_to_plot:
                if key == par_name:
                    names.append(key)
                    labels.append(par['label'])
                    if 'log' in par:
                        samps[:,j] = np.log10(param_values[key])
                        ranges.update({key:np.log10(par['range'])})
                    else:
                        samps[:,j] = param_values[key]
                        ranges.update({key:par['range']})
                    j += 1
                    break

        settings = {
                    "contours":[0.68, 0.95, 0.99], 
                    "range_ND_contour":1, 
                    "range_confidence":0.001,
                    "fine_bins":200,
                    "fine_bins_2d":80,
                    "smooth_scale_1D":0.3,
                    "smooth_scale_2D":0.5,
                    "tight_gap_fraction":0.15
                    }

        samples = MCSamples(samples=samps, names=names, ranges=ranges, \
                weights=probability, labels=labels, settings=settings )

        g = plots.getSubplotPlotter()
        g.settings.num_plot_contours = 3
        g.settings.prob_y_ticks = True

        line_args = {"lw":2, "color":colorConverter.to_rgb("#006FED") } 

        g.triangle_plot(samples, filled=True, line_args=line_args)

        g.fig.subplots_adjust(wspace=0.1, hspace=0.1)

        # Add tick labels at top of diagonal panels
        for i, ax in enumerate([g.subplots[i,i] for i in range(nParamsToPlot)]):
            par_name = keys[i]
            if i < nParamsToPlot-1: 
                ticklabels = [item.get_text() for item in g.subplots[-1,i].xaxis.get_ticklabels()]
                ax.xaxis.set_ticklabels(ticklabels)
                ax.tick_params(labelbottom='off',labeltop='on')

            # Add shaded region showing 1D 68% credible interval
            y0, y1 = ax.get_ylim()
            lev = samples.get1DDensity(par_name).getLimits(settings['contours'][0])

            ax.add_patch(
                    Rectangle((lev[0], y0), 
                    lev[1]-lev[0], 
                    y1-y0, 
                    facecolor="grey", 
                    alpha=0.5)
                    )

            # Indicate the value of the "true" parameter
            if mock_param_values is not None:
                ax.plot(
                        mock_param_values[i],
                        y0+(y1-y0)*0.05,
                        marker="D",
                        ms=8,
                        color="green") 

        if show:
            plt.show()
        else:
            # Now save the plot
            name = prepare_plot_saving(plot_name)
            g.export(name)

        plt.close()
        hdulist.close()


##                # Overplot the posterior median point
##                self.marginal.DrawPosteriorMedian(Appearance.PosteriorMedian)
##
##                # Overplot the posterior mean point
##                self.marginal.DrawPosteriorMean(Appearance.PosteriorMean)
##
##                # Overplot the maximum a posteriori
##                self.marginal.DrawMAP(Appearance.BestFit)
##
##                # Overplot the 1-D credible regions
##                self.marginal.DrawCredibleRegions(
##                    self.CredibleRegions)
##
##                # If has mock parameters, overplot the "true" values
##                try:
##                    mockValue = self.mockpar[iRow]
##                    if self.logPar[iRow]:
##                        mockValue = np.log10(mockValue)
##                    self.OneDim.PlotMockPar(
##                            mockValue,
##                            Appearance.MockPar)
##                except AttributeError:
##                    None
##
##                # Overplot the posterior mean point
##                self.joint.DrawPosteriorMean(
##                    Appearance.PosteriorMean)

        #plt.close(g.fig)
