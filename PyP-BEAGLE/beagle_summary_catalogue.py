import os
import logging
from collections import OrderedDict
import json
import numpy as np
from scipy.integrate import cumtrapz
from scipy.interpolate import interp1d
from astropy.io import fits

from beagle_utils import prepare_data_saving, BeagleDirectories
from significant_digits import to_precision, float_nsf

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

    # ******************************************************************
    # Here you must simply use `cumsum`, and not `cumtrapz` as in
    # beagle_utils.prepare_violin_plot, since the output of MultiNest are a set
    # of weights (which sum up to 1) associated to each set of parameters (the
    # `p_j` of equation 9 of Feroz+2009), and not a probability density (as the
    # MultiNest README would suggest).
    # ******************************************************************
    cumul_pdf = np.cumsum(probability[sort_])
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


class BeagleSummaryCatalogue(object):


    def load(self, file_name=None):
        """ 
        Load a 'BEAGLE summary catalogue'

        Parameters
        ----------
        file_name : str
            Name of the file containing the catalogue.
        """

        if file_name is None:
            file_name = "BEAGLE_summary_catalogue.fits"

        logging.info("Loading the `BeagleSummaryCatalogue` file: " + file_name)

        name = file_name
        if not os.path.dirname(file_name):
            name = os.path.join(BeagleDirectories.results_dir, 
                    BeagleDirectories.pypbeagle_data, 
                    file_name)

        self.hdulist = fits.open(name)

    def compute(self, file_list, config_file, file_name=None, overwrite=False, levels=[68.,95.]):
        """ 
        """ 

        if file_name is None:
            file_name = "BEAGLE_summary_catalogue.fits"

        with open(config_file) as f:    
            hdu_col = json.load(f, object_pairs_hook=OrderedDict)


        # You consider the first file in the list and use as a "mold" to create
        # the structure (binary tables and their columns) of the output FITS file
        firstfile = os.path.join(BeagleDirectories.results_dir, file_list[0])
        hdulist = fits.open(firstfile)

        n_objects = len(file_list)

        # Initialize a new (empty) primary HDU for your output FITS file
        self.hdulist = fits.HDUList(fits.PrimaryHDU())
    
        # Now you cycle over all extension and columns that you want to put in
        # the summary catalogue
        for hdu in hdu_col:
            new_columns = list()

            # Pick the extension name
            hdu_name = hdu['name']

            # The first column of each output extension contains the object ID
            #new_columns.append(fits.Column(name='ID', format='K'))
            new_columns.append(fits.Column(name='ID', format='20A'))

            # You just consider the columns defined in the structure
            if 'columns' in hdu:
                columnNames = hdu['columns']
            # While by default you take all columns in that extensions
            else:
                columnNames = hdulist[hdu_name].columns.names

            # For each column, you add a '_mean', '_median' and confidence
            # intervals columns, taking the appropriate units from the FITS
            # file that you are using as a mold
            for col_index in range(len(hdulist[hdu_name].columns)):

                col_ = hdulist[hdu_name].columns[col_index]

                if col_.name not in columnNames:
                    continue
    
                new_columns.append(fits.Column(name=col_.name+'_mean',
                    format=col_.format, unit=col_.unit))

                new_columns.append(fits.Column(name=col_.name+'_median',
                    format=col_.format, unit=col_.unit))

                for lev in np.array(levels):
                    new_columns.append(fits.Column(name=col_.name + '_' +
                        "{:.2f}".format(lev), format='2'+col_.format[-1],
                        unit=col_.unit))

            # Create the "column definition"
            cols_ = fits.ColDefs(new_columns)

            # Create the actual binary table, with the correct number of rows
            # to accomodate all objects
            nrows = len(file_list)
            new_hdu = fits.BinTableHDU.from_columns(cols_, nrows=nrows)
            new_hdu.name = hdu_name

            # And finally append the newly created bunary table to the hdulist
            # that will be printed to the ouput FITS file
            self.hdulist.append(new_hdu)

        hdulist.close()

        # Now you can go through each file, and compute the required quantities

        for i, file in enumerate(file_list):
            hdulist = fits.open(os.path.join(BeagleDirectories.results_dir, file))
            end = file.find('_BEAGLE')

            # Extract the object ID from the file_name
            #ID = np.int(np.float(os.path.basename(file[0:end])))
            ID = os.path.basename(file[0:end])

            probability = hdulist['posterior pdf'].data['probability']

            for hdu in hdu_col:
                hdu_name = hdu['name']
                if 'columns' in hdu:
                    columnNames = hdu['columns']
                else:
                    columnNames = hdulist[hdu_name].columns.names

                for col_name in columnNames:
                    self.hdulist[hdu_name].data['ID'][i] = ID
                    par_values = hdulist[hdu_name].data[col_name]

                    mean, median, interval = get1DInterval(par_values, probability, levels)

                    self.hdulist[hdu_name].data[col_name+'_mean'][i] = mean
                    self.hdulist[hdu_name].data[col_name+'_median'][i] = median

                    for j, lev in enumerate(levels):
                        levName = col_name + '_' + "{:.2f}".format(lev)
                        self.hdulist[hdu_name].data[levName][i] = interval[j]

            hdulist.close()

        name = prepare_data_saving(file_name)
        self.hdulist.writeto(name, clobber=overwrite)



    def make_latex_table(self, param_names, 
            IDs=None,
            level=68., 
            summary_statistics='median',
            significant_digits=2):

        n = significant_digits
        
        n_rows = len(self.hdulist[1].data['ID'])

        if IDs is None:
            IDs = self.hdulist[1].data['ID']

        for ID in IDs:
            print "\nID: ", ID, 
            row = np.arange(n_rows)[self.hdulist[1].data['ID'] == ID][0]
            for param in param_names:
                colName = param + '_' + summary_statistics
                for hdu in self.hdulist:
                    if hasattr(hdu, 'columns'):
                        if colName in hdu.columns.names:
                            sum_stat = hdu.data[colName][row]
                            errColName = param + '_' + "{:.2f}".format(level)
                            error = hdu.data[errColName][row]

                            print to_precision(sum_stat, n+2), '_{-' + to_precision(sum_stat-error[0], n) \
                                    + '}^{+' + to_precision(error[1]-sum_stat, n) + '} &',
                            break

                        
