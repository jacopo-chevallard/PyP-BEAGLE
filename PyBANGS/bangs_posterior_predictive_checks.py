import os
import numpy as np
from astropy.io import fits
from astropy.table import Table, Column

from bangs_utils import prepare_file_writing

class PosteriorPredictiveChecks:


    def load(self, file_name):
        """ 
        Load a file containing the posterior predictive checks

        Parameters
        ----------
        file_name : str
            Name of the file.
        """

        my_table = Table.read(file_name)
    
        self.data = my_table

    def compute(self, observed_catalogue, filters, 
            results_dir, file_name=None):
        """ 
        Compute  posterior predictive checks quantities.

        Parameters
        ----------
        observed_catalogue : `bangs_photometry.ObservedCatalogue`
            Class containing an observed photometric catalogue.

        filters : `bangs_filters.PhotometricFilters`
            Class containing a set of photometric filters.

        results_dir : str
            Directory containing the BANGS output files.

        file_name : str, optional
            Name of the output catalogue, wuthout including the direcory tree.
            It will be saved into the `results_dir`/pybangs/data folder (which
            will be created if not present).
        """

##        if filters is not None:
##            self.filters = filters
##
##        if not hasattr(self, 'filters'):
##            raise AttributeError("No 'PhotometricFilters' class has been passed"
##                    " to the class constructor or to this function!")
##
##        if results_dir is not None:
##            self.results_dir = results_dir
##
##        if not hasattr(self, 'results_dir'):
##            raise AttributeError("No 'results_dir' defined!")

        # Copy from the catalogue the column containing the object IDs
        objID = Column(data=observed_catalogue.data['ID'], name='ID', dtype=np.int32) 

        n_obj = len(observed_catalogue.data['ID'])
        
        # Defines columns containing the number of photometric bands actually
        # used in the BANGS run, for a given object,
        n_used_bands = Column(name='n_used_bands', dtype=np.int32, length=n_obj)

        # the average chi-square, and average chi-square/n_used_bands
        aver_chi_square = Column(name='aver_chi_square', dtype=np.float32,
                length=n_obj)

        aver_red_chi_square = Column(name='aver_red_chi_square',
                dtype=np.float32, length=n_obj)

        my_cols = [objID, n_used_bands, aver_chi_square, aver_red_chi_square]
        my_table = Table(my_cols)

        obs_flux = np.zeros(filters.n_bands, np.float32)
        obs_flux_err = np.zeros(filters.n_bands, np.float32)

        model_flux = np.zeros(filters.n_bands, np.float32)

        jy = 1.E-26

        for i in range(n_obj):

            strID = str(objID[i])
            file = os.path.join(results_dir, strID + "_BANGS.fits.gz")

            if os.path.isfile(file):

                # Open the FITS file containing BANGS results for the current object
                hdulist = fits.open(file)
                bangs_data = hdulist['MARGINAL PHOTOMETRY'].data

                probability = hdulist['POSTERIOR PDF'].data['probability']

                n_samples = len(bangs_data.field(0))
                chi_square = np.zeros(n_samples, np.float32)
                n_data = 0

                for j in range(filters.n_bands):

                    # observed flux and its error
                    name = filters.data['flux_colName'][j]
                    obs_flux = observed_catalogue.data[i][name] * filters.units / jy

                    name = filters.data['flux_errcolName'][j]
                    obs_flux_err = observed_catalogue.data[i][name] * filters.units  / jy

                    # if defined, add the minimum error in quadrature
                    obs_flux_err = np.sqrt((obs_flux_err/obs_flux)**2 + np.float32(filters.data['min_rel_err'][j])**2) * obs_flux

                    # model flux and its error
                    name = '_' + filters.data['label'][j] + '_'
                    model_flux = bangs_data[name] / jy

                    if obs_flux_err > 0.:
                        n_data += 1
                        chi_square += ((obs_flux-model_flux) / obs_flux_err)**2

                my_table['n_used_bands'][i] = n_data
                my_table['aver_chi_square'][i] = np.sum(probability*chi_square) / np.sum(probability)
                my_table['aver_red_chi_square'][i] = my_table['aver_chi_square'][i] / (n_data-1)

                hdulist.close()

        self.columns = my_cols
        self.data = my_table

        if file_name is not None:
            name = prepare_file_writing(results_dir, file_name)
            my_table.write(name)

