import numpy as np
import os
from astropy.table import Table, Column
from astropy.io import fits

class PosteriorPredictiveChecks:

    def load_catalogue(self, tablename):

        try:
            my_table = Table.read(tablename)
        except IOError as e:
            print "I/O error({0}): {1}".format(e.errno, e.strerror)
    
        self.data = my_table

    def make_catalogue(self, observed_catalogue, resultsDir, filters, tableName=None):

        objID = Column(data=observed_catalogue['ID'], name='ID', dtype=np.int32) 

        n_obj = len(observed_catalogue['ID'])

        n_used_bands = Column(name='n_used_bands', dtype=np.int32, length=n_obj)
        aver_chi_square = Column(name='aver_chi_square', dtype=np.float32, length=n_obj)
        aver_red_chi_square = Column(name='aver_red_chi_square', dtype=np.float32, length=n_obj)

        my_cols = [objID, n_used_bands, aver_chi_square, aver_red_chi_square]
        my_table = Table(my_cols)

        obs_flux = np.zeros(filters.n_bands, np.float32)
        obs_flux_err = np.zeros(filters.n_bands, np.float32)

        model_flux = np.zeros(filters.n_bands, np.float32)

        jy = 1.E-26

        for i in range(n_obj):

            strID = str(objID[i])
            file = os.path.join(resultsDir, strID + "_BANGS.fits.gz")

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
                    obs_flux = observed_catalogue[i][filters.colName[j]] * filters.units / jy
                    obs_flux_err = observed_catalogue[i][filters.errcolName[j]] * filters.units  / jy

                    # if defined, add the minimum error in quadrature
                    if filters.min_rel_err:
                        obs_flux_err = np.sqrt((obs_flux_err/obs_flux)**2 + np.float32(filters.min_rel_err[j])**2) * obs_flux

                    # model flux and its error
                    name = '_' + filters.label[j] + '_'
                    model_flux = bangs_data[name] / jy

                    if obs_flux_err > 0.:
                        n_data += 1
                        chi_square += ((obs_flux-model_flux) / obs_flux_err)**2

                my_table['n_used_bands'][i] = n_data
                my_table['aver_chi_square'][i] = np.sum(probability*chi_square) / np.sum(probability)
                my_table['aver_red_chi_square'][i] = my_table['aver_chi_square'][i] / n_data

                hdulist.close()

        self.data = my_table

        if tableName is not None:
            directory = os.path.join(resultsDir, 'pybangs', 'data')
            if not os.path.exists(directory):
                os.makedirs(directory)

            name = os.path.join(directory, tableName)
            my_table.write(name)

