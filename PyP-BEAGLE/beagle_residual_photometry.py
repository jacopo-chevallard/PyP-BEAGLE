from __future__ import absolute_import
from __future__ import print_function
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import six.moves.cPickle

import os
import sys
import pyp_beagle.dependencies.WeightedKDE 

from .beagle_utils import weighted_avg_and_std, BeagleDirectories, match_ID, prepare_plot_saving
from six.moves import range

class ResidualPhotometry(object):

    def load(self, file_name):
        """ 
        Load the pre-computed residul photometry from a cPickle file.

        Parameters
        ----------
        file_name : str
            Name of the file containing the cPicke dump.
        """ 

        name = os.path.join(BeagleDirectories.results_dir,
                BeagleDirectories.pypbeagle_data, file_name)

        file = open(name, 'rb')
        self.residual_kde_pdf = six.moves.cPickle.load(file)
        file.close()

    def compute(self, observed_catalogue, beagle_summary_catalogue,
            filters, summary_stat=None, cPickleName=None):

        """ 
        Compute the residual photometry between an observed and model catalogues.  

        Parameters
        ----------
        observed_catalogue : `ObservedCatalogue` class object
            

        beagle_summary_catalogue : `BeagleSummaryCatalogue` class object
            

        filters : `PhotometricFilters` class object
            
        summary_stat : str, optional
            Either 'mean' or 'median', determines which type of summary
            statistics is used to compute the residual. By default the median
            of the marginal PDF is used.

        """

        if summary_stat is None:
            summary_stat = "median"

        beagle_data = beagle_summary_catalogue.hdulist['MARGINAL PHOTOMETRY'].data

        catalogue_data = observed_catalogue.data

        indx_beagle, indx_catalogue = match_ID(beagle_data['ID'], catalogue_data['ID'])

        beagle_data = beagle_data[indx_beagle]
        catalogue_data = catalogue_data[indx_catalogue]

        print("ID: ", beagle_data['ID'])
        print("ID: ", catalogue_data['ID'])

        # As a sanity check, check if the ID match among the two catalogues, by
        # random picking some indices here and there...
        #if (beagle_data['ID'][0] != catalogue_data['ID'][0]) or \
        #(beagle_data['ID'][-1] != catalogue_data['ID'][-1]):
        #    raise ValueError("The object IDs between the BEAGLE summary catalogue and \
        #        the observed catalogue do not match!")

        self.residual_kde_pdf = list() 

        jy = 1.E-26

        for i in range(filters.n_bands):

            obs_flux = catalogue_data[filters.data['flux_colName'][i]] * filters.units / jy
            obs_flux_err = catalogue_data[filters.data['flux_errcolName'][i]] * filters.units  / jy

            name = '_' + filters.data['label'][i] + '_'
            model_flux = beagle_data[name+'_'+summary_stat] / jy
            model_flux_err = 0.5 * (beagle_data[name+'_68.00'][:,1]-beagle_data[name+'_68.00'][:,0]) / jy

            mask = np.zeros(len(obs_flux), dtype=bool)

            mask[(obs_flux > 0.) & (obs_flux_err > 0.) & (model_flux > 0.) & (model_flux_err > 0.)] = True

            # Compute the residual in magnitude scale
            residual = -2.5*np.log10(model_flux[mask]/obs_flux[mask])

            # The residual error follows from standard error propagation
            residualErr = np.sqrt((2.5/model_flux[mask]*model_flux_err[mask])**2 + (2.5/obs_flux[mask]*obs_flux_err[mask])**2)

            # Compute the distribution of the residual by means of a weighted kernel density estimation
            # Sometimes at high redshift, short wavelength filters all contain model_flux = 0 when wavelength
            # of filter outside of wavelength range of templates, hence testing the length of residual
            if len(residual) > 0:
                kde_pdf = WeightedKDE.gaussian_kde(residual, weights = 1./residualErr)
            else:
                kde_pdf = 0

            self.residual_kde_pdf.append(kde_pdf)

        if cPickleName is not None:
            six.moves.cPickle.dump(self.residual_kde_pdf, cPickleName, six.moves.cPickle.HIGHEST_PROTOCOL)

    def plot(self, 
            filters, 
            plot_name="BEAGLE_residual_photometry.pdf", 
            x_range=(-1.,1.), 
            n_x=1000, 
            summary_stat="median"):
        """ 
        Plot the residual photometry.

        Parameters
        ----------
        filters : PhotometricFilters class object

        plot_name: str, optional
            File name of the output plot.

        x_range : iterable float size=2, optional
            Range over which computing the actual residual density function
            from the computed KDE residual.

        n_x : int, optional 
            Number of values over which computing the residual density
            function.

        summary_stat : str, optional
            Either 'mean' or 'median', determines which type of summary
            statistics is used.
            
        Notes
        -----
        The residuals density function is computed at each residual value
        defined by np.linspace(x_range[0], x_range[1], n_x)
        """

        # The grid of residual values over which computing the residual density
        # function
        x_grid = np.linspace(x_range[0], x_range[1], n_x)

        for residual_kde_pdf in self.residual_kde_pdf:

            # The residual density function, computed from the residual KDE
            kde_pdf_grid = np.array(residual_kde_pdf(x_grid))

            # You also compute the median, or mean residual, and its dispersion (or
            # 68 % credible region)
            cumul_pdf = np.cumsum(kde_pdf_grid)

            # Be sure the it sums up to 1
            cumul_pdf /= cumul_pdf[-1]

            # Compute interpolant of cumulative PDF (linear interpolation)
            interp_cumul_pdf = interp1d(cumul_pdf, x_grid)

            if "median" in summary_stat:
                med = interp_cumul_pdf(0.5)
                # This corresponds to a 68 % credible region
                interval = interp_cumul_pdf(0.84) - interp_cumul_pdf(0.16)
                print("\n median = {:.3f}".format(med))
                print("68 % interval = {:.3f}".format(interval))
            elif "mean" in summary_stat:
                mean, stddev = weighted_avg_and_std(x_grid, kde_pdf_grid)
                print("\n mean = {:.3f}".format(mean))
                print("stddev = {:.3f}".format(stddev))

        # Compute interpolant of cumulative PDF (linear interpolation)
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)

        ax.plot(x_grid, kde_pdf_grid , ls="-", color = "black")

        name = prepare_plot_saving(plot_name)

        fig.savefig(name, dpi=None, facecolor='w', edgecolor='w',
            orientation='portrait', papertype='a4', format="pdf",
            transparent=False, bbox_inches="tight", pad_inches=0.1)

        plt.close(fig)


