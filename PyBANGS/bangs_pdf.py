import json
from astropy.io import fits

from getdist import plots, MCSamples

class PDF:

    def __init__(self, param_names_file):

        # Names of parameters, used to label the axes, and whether they are log
        # or not
        with open(param_names_file) as f:    
            self.adjust_params = json.load(f)

    def plot_triangle(self, fits_file, params_to_plot=None, suffix=None):

        hdulist = fits.open(fits_file)

        param_values = hdulist['posterior pdf'].data
        probability = hdulist['posterior pdf'].data['probability']

        n_rows = probability.size

       # ParamsToPlot = ['mass', 'redshift', 'tauV_eff', 'metallicity', 'sfr', 'tau']

        # By default you plot all parameters
        if params_to_plot is None:
            _params_to_plot = 
        else: 
            _params_to_plot = params_to_plot

        nParamsToPlot = len(ParamsToPlot)

        names = list()
        labels = list()
        ranges = dict()
        samps = np.zeros((n_rows, nParamsToPlot))

        j = 0
        for par_name in ParamsToPlot:
            for key, par in self.adjust_params.iteritems():
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
                    "smooth_scale_1D":0.5,
                    "smooth_scale_2D":0.7
                    }

        samples = MCSamples(samples=samps, names=names, ranges=ranges,
                weights=probability, labels=labels, settings=settings )

        g = plots.getSubplotPlotter()
        g.settings.num_plot_contours = 3
        g.settings.prob_y_ticks = True

        line_args = {"lw":2, "color":colorConverter.to_rgb("#006FED") } 

        g.triangle_plot(samples, filled=True, line_args=line_args)

        g.fig.subplots_adjust(wspace=0.1, hspace=0.1)

        # Add tick labels at top of diagonal panels
        for i, ax in enumerate([g.subplots[i,i] for i in range(nParamsToPlot)]):
            par_name = ParamsToPlot[i]
            if i < nParamsToPlot-1: 
                ticklabels = [item.get_text() for item in g.subplots[-1,i].xaxis.get_ticklabels()]
                ax.xaxis.set_ticklabels(ticklabels)
                ax.tick_params(labelbottom='off',labeltop='on')

            # Add shaded region showing 1D 68% credible interval
            y0, y1 = ax.get_ylim()
            lev = samples.get1DDensity(par_name).getLimits(settings['contours'][0])
            ax.add_patch(Rectangle((lev[0], y0), lev[1]-lev[0], y1-y0, facecolor="grey", alpha=0.2))

        # Now save the plot
        if suffix is None:
            outputName = self.outputFileName + '_triangle.pdf'
        else:
            outputName = self.outputFileName + '_triangle_' + suffix + '.pdf'

        g.export( outputName )

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

        plt.close(g.fig)
