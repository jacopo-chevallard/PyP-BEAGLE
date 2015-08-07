import numpy as np
import os
from astropy.io import fits
from scipy.integrate import cumtrapz
from scipy.interpolate import interp1d


def get1DInterval(param_values, probability, levels):
        """ 
        Compute several quantities from a 1D probability density function

        Parameters
        ----------
        param_values : numpy array
            Contains the values of the parameter.

        probability : numpy array 
            Contains the probability associated with each value of the
            parameter, hence must have same dimension as 'param_values'.

        levels : numpy array or list containing float
            Contains the (percentage) levels used to compute the credible
            regions, e.g. levels=[68.,95.] will compute 68 % and 95 % (central)
            credible regions

        Returns
        -------
        mean : float
            Mean of the parameter, computed as
            sum(probability*param_values)/sum(probability)
        median : float
            Median of the parameter, computed from the cumulative integral of
            the PDF
        interval : list of float
            2-dimensional list containing the lower and upper value of the
            parameter corresponding to the different `levels`

        """

    sort_ = np.argsort(param_values)
    cumul_pdf = cumtrapz(probability[sort_], param_values[sort_], initial = 0.)
    cumul_pdf /= cumul_pdf[len(cumul_pdf)-1]

    # Get the interpolant of the cumulative probability
    f_interp = interp1d(cumul_pdf, param_values[sort_])

    # You shoud integrate rather than summing here
    mean = np.sum(probability * param_values) / np.sum(probability)

    median = f_interp(0.5)

    interval = list()
    for lev in levels:

        low, high = f_interp([0.5*(1.-lev/100.), 1.-0.5*(1.-lev/100.)])
        interval.append([low,high])

    return mean, median, interval


class BangsSummaryCatalogue:


    def load_catalogue(self, filename):

        try:
            self.hdulist = fits.open(filename)
        except IOError as e:
            print "I/O error({0}): {1}".format(e.errno, e.strerror)


    def make_catalogue(self, filelist, filename, levels=[68.,95.]):

        hdu_col = list()

        hdu_col.append({'name':'GALAXY PROPERTIES', 'columns':['redshift',
            'M_tot', 'M_star', 'mass_w_age', 'lumin_w_age', 'mass_w_Z',
            'lumin_w_Z', 'N_ion', 'xi_ion', 'UV_slope']})

        hdu_col.append({'name':'STAR FORMATION', 'columns':['SFR', 'sSFR']})

        hdu_col.append({'name':'DUST ATTENUATION', 'columns':['tauV_eff']})

        hdu_col.append({'name':'ABSOLUTE MAGNITUDES'})

        # You considerthe first file in the list and use as a "mold" to create
        # the structure (binary tables and their columns) of the output FITS file
        firstfile = filelist[0]
        hdulist = fits.open(firstfile)

        n_objects = len(filelist)

        # Initialize a new (empty) primary HDU for your output FITS file
        my_hdu = fits.HDUList(fits.PrimaryHDU())
    
        # Now you cycle over all extension and columns that you want to put in
        # the summary catalogue
        for hdu in hdu_col:
            new_columns = list()

            # Pick the extension name
            hduName = hdu['name']

            # The first column of each output extension contains the object ID
            new_columns.append(fits.Column(name='ID', format='20A'))

            # You just consider the columns defined in the structure
            if 'columns' in hdu:
                columnNames = hdu['columns']
            # While by default you take all columns in that extensions
            else:
                columnNames = hdulist[hduName].columns.names

            # For each column, you add a '_mean', '_median' and confidence
            # intervals columns, taking the appropriate units from the FITS
            # file that you are using as a mold
            for colName in columnNames:
                col_ = hdulist[hduName].columns[colName]

                new_columns.append(fits.Column(name=col_.name+'_mean',
                    format=col_.format, unit=col_.unit))

                new_columns.append(fits.Column(name=col_.name+'_median',
                    format=col_.format, unit=col_.unit))

                for lev in levels:
                    new_columns.append(fits.Column(name=col_.name + '_' +
                        "{:.2f}".format(lev), format='2'+col_.format[-1],
                        unit=col_.unit))

            # Create the "column definition"
            cols_ = fits.ColDefs(new_columns)

            # Create the actual binary table, with the correct number of rows
            # to accomodate all objects
            new_hdu = fits.BinTableHDU.from_columns(cols_, nrows=nrows)
            new_hdu.name = hduName

            # And finally append the newly created bunary table to the hdulist
            # that will be printed to the ouput FITS file
            my_hdu.append(new_hdu)

        hdulist.close()

        # Now you can go through each file, and compute the required quantities

        for i, file in enumerate(filelist):
            hdulist = fits.open(file)
            end = BANGS_file.find('_BANGS')

            # Extract the object ID from the filename
            ID = os.path.basename(BANGS_file[0:end])

            probability = hdulist['posterior pdf'].data['probability']

            for hdu in hdu_col:
                hduName = hdu['name']
                if 'columns' in hdu:
                    columnNames = hdu['columns']
                else:
                    columnNames = hdulist[hduName].columns.names

                for colName in columnNames:
                    my_hdu[hduName].data['ID'][i] = ID
                    par_values = hdulist[hduName].data[colName]

                    mean, median, interval = get1DInterval(par_values, probability, levels)

                    my_hdu[hduName].data[colName+'_mean'][i] = mean
                    my_hdu[hduName].data[colName+'_median'][i] = median

                    for j, lev in enumerate(levels):
                        levName = colName + '_' + "{:.2f}".format(lev)
                        my_hdu[hduName].data[levName][i] = interval[j]

            hdulist.close()

        my_hdu.writeto(filename, clobber=True)
