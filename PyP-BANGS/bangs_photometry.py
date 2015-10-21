import os
from scipy.integrate import simps, cumtrapz
from scipy.interpolate import interp1d
from bisect import bisect_left
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as plticker
from astropy.io import ascii
from astropy.io import fits

import sys
sys.path.append("../dependencies")
import WeightedKDE

from bangs_utils import BangsDirectories, prepare_plot_saving, set_plot_ticks
from bangs_filters import PhotometricFilters
from bangs_summary_catalogue import BangsSummaryCatalogue
from bangs_residual_photometry import ResidualPhotometry
from bangs_multinest_catalogue import MultiNestCatalogue
from bangs_posterior_predictive_checks import PosteriorPredictiveChecks


microJy = np.float32(1.E-23 * 1.E-06)
nanoJy = np.float32(1.E-23 * 1.E-09)

class ObservedCatalogue:

    def load(self, file_name):

        """ 
        Load a photometric catalogue of observed sources. It automatically
        detects, and loads, FITS or ASCII files depending on the suffix.

        Parameters
        ----------
        file_name : str
            Contains the file name of the catalogue.
        """

        if file_name.endswith(('fits', 'fit', 'FITS', 'FIT')):
            self.data = fits.open(file_name)[1].data
            self.columns = fits.open(file_name)[1].columns
        else:
            self.data = ascii.read(file_name, Reader=ascii.basic.CommentedHeader)


class Photometry:

    def __init__(self):

        self.filters = PhotometricFilters()

        self.observed_catalogue = ObservedCatalogue()

        self.summary_catalogue = BangsSummaryCatalogue()

        self.multinest_catalogue = MultiNestCatalogue()

        self.residual = ResidualPhotometry()

        self.PPC = PosteriorPredictiveChecks()

    def plot_marginal(self, ID, max_interval=99.7, print_text=False):    
        """ 
        Plot the 'marginal' photometry predicted by BANGS.

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
        """

        # From the (previously loaded) observed catalogue select the row
        # corresponding to the input ID
        observation = self.observed_catalogue.data[self.observed_catalogue.data['ID'] == ID]

        # Check if you need to apply an aperture correction to the catalogue fluxes
        if 'aper_corr' in self.observed_catalogue.data.dtype.names:
            aper_corr = 10.**(-0.4*observation[0]['aper_corr'])
        else:
            aper_corr = 1.

        # Put observed photometry and its error in arrays
        obs_flux = np.zeros(self.filters.n_bands)
        obs_flux_err = np.zeros(self.filters.n_bands)

        for i, band in enumerate(self.filters.data['flux_colName']):
            obs_flux[i] = observation[0][band]*aper_corr*self.filters.units / nanoJy

        # Add to the error array the minimum relative error thet BANGS allows
        # one to add to the errors quoted in the catalogue
        for i, err in enumerate(self.filters.data['flux_errcolName']):
            tmp_err = observation[0][err]
            if tmp_err > 0.:
                obs_flux_err[i] = observation[0][err]*aper_corr*self.filters.units / nanoJy
                obs_flux_err[i] = (np.sqrt( (obs_flux_err[i]/obs_flux[i])**2 +
                        np.float32(self.filters.data['min_rel_err'][i])**2) *
                        obs_flux[i])
            else:
                obs_flux_err[i] = tmp_err

        ok = np.where(obs_flux_err > 0.)[0]

        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)

        # Open the file containing BANGS results
        fits_file = os.path.join(BangsDirectories.results_dir, str(ID)+'_BANGS.fits.gz')
        hdulist = fits.open(fits_file)

        # Consider only the extension containing the predicted model fluxes
        model_sed = hdulist['marginal photometry']
        probability = hdulist['posterior pdf'].data['probability']

        n_bands = len(model_sed.columns.names)
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

        width = 5*np.min(np.array(self.filters.data['wl_eff'][1:])-np.array(self.filters.data['wl_eff'][0:-1]))

        kwargs = {'color':'tomato', 'alpha':0.7, 'edgecolor':'black', 'linewidth':0.2}

        for i in range(n_bands):
            xdata = model_sed.data.field(i) / nanoJy

            min_x = np.min(xdata)
            max_x = np.max(xdata)

            min_flux[i] = min_x
            max_flux[i] = max_x

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
            pdf = WeightedKDE.gaussian_kde(xdata, weights=probability)

            # Build a grid of value over which computing the actual PDF from its KDE
            x_grid = np.linspace(min_x, max_x, nXgrid)

            # Compute the PDF
            pdf_grid = np.array(pdf(x_grid))

            # Compute the PDF normalization, so its integral will be exactly 1
            pdf_norm[i] = simps(pdf_grid, x_grid)
            pdf_grid /= pdf_norm[i]

            # Now compute the cumulative PDF
            cumul_pdf = cumtrapz(pdf_grid, x_grid, initial = 0.)
            cumul_pdf /= cumul_pdf[-1]

            # Get the interpolant of the cumulative PDF
            f_interp = interp1d(cumul_pdf, x_grid)

            # Compute the limits over which you will plot the "violin", for
            # instance showing only the cumulative PDF up to +/- 3 sigma
            intv = 0.5*(100.-max_interval)/100.
            lims = f_interp([intv,1.-intv])
            prob_lims = pdf(lims) / pdf_norm[i]

            # The median corresponds to a cumulative probability = 0.5
            median_flux[i] = f_interp(0.5)

            kde_pdf[i] = pdf

            i1 = bisect_left(x_grid, lims[0])
            i2 = bisect_left(x_grid, lims[1])
        
            x_plot[i] = np.concatenate(([lims[0]], x_grid[i1+1:i2], [lims[1]]))

            y_plot[i] = np.concatenate(([prob_lims[0]], pdf_grid[i1+1:i2], [prob_lims[1]]))

            _max_y[i] = np.max(y_plot[i])

        delta_wl = np.array(self.filters.data['wl_eff'][1:])-np.array(self.filters.data['wl_eff'][0:-1])
        delta_wl = np.concatenate(([delta_wl[0]], delta_wl))

        for i in range(n_bands):

            dwl = delta_wl[i]
            if i > 1:
                dwl = np.min(delta_wl[i-1:i])

            if has_pdf[i]:

                w = 0.4 * dwl / _max_y[i]

                y_grid = np.full(len(x_plot[i]), self.filters.data['wl_eff'][i])
                _lim_y = kde_pdf[i](median_flux[i])/pdf_norm[i] * w

                ax.fill_betweenx(x_plot[i],
                        y_grid - y_plot[i]*w,
                        y_grid + y_plot[i]*w,
                        **kwargs
                        )

                ax.plot( [self.filters.data['wl_eff'][i]-_lim_y, self.filters.data['wl_eff'][i]+_lim_y],
                        [median_flux[i], median_flux[i]],
                        color = 'black',
                        linewidth = 0.2
                        )

            ax.plot( self.filters.data['wl_eff'][i],
                    median_flux[i],
                    color = 'black',
                    marker = "o",
                    markersize = 5,
                    alpha = 0.7
                    )


        # Determine min and max values of y-axis
        yMax = np.max(max_flux)
        yMin = np.min(np.concatenate((obs_flux[ok]-obs_flux_err[ok], min_flux)))

        dY = yMax-yMin

        yMax += dY * 0.1
        yMin -= dY * 0.1

        ax.set_ylim([yMin, yMax])

        x0 = self.filters.data['wl_eff'][0]
        x1 = self.filters.data['wl_eff'][-1]
        dx = x1-x0
        ax.set_xlim([x0-0.05*dx, x1+0.05*dx])

        x0, x1 = ax.get_xlim()
        if yMin < 0.: plt.plot( [x0,x1], [0.,0.], color='gray', lw=0.8 )

        # Define plotting styles
        ax.set_xlabel("$\lambda_\\textnormal{eff} / \\textnormal{\AA}$ (observed-frame)")
        ax.set_ylabel("$f_{\\nu}/\\textnormal{nanoJy}$")

        # Set better location of tick marks
        set_plot_ticks(ax, n_x=5)

        kwargs = {'alpha':0.7}

        plt.errorbar(self.filters.data['wl_eff'], 
                obs_flux, 
                yerr = obs_flux_err,
                color = "dodgerblue",
                ls = " ",
                marker = "D",
                markeredgewidth = 0.,
                markersize = 8,
                elinewidth=0.5,
                **kwargs)


        # Title of the plot is the object ID
        plt.title(str(ID))

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

        name = prepare_plot_saving(str(ID)+'_BANGS_marginal_SED_phot.pdf')

        fig.savefig(name, dpi=None, facecolor='w', edgecolor='w',
                orientation='portrait', papertype='a4', format="pdf",
                transparent=False, bbox_inches="tight", pad_inches=0.1)
        plt.close(fig)

        hdulist.close()

##    def plot_residuals(self, residual_file_name=None, residual_plotname=None):
##
##        if not hasattr(self, 'observed_catalogue'):
##            except AttributeError:
##                    "An observed catalogue must be loaded before plotting the
##                    residual"
##
##        if not hasattr(self, 'bangs_summary_catalogue'):
##            except AttributeError:
##                    "A `bangs_summary_catalogue` must be loaded before plotting the
##                    residual"
##
##        self.residual = ResidualPhotometry()
##
##        try:
##            self.residual.load(self.residual_file_name)
##        except:
##            self.residual.compute(self.observed_catalogue,
##                self.bangs_summary_catalogue, self.self.filters.
##                cPickleName=self.residual_file_name)
                


