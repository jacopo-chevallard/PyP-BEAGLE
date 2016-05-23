import os
import logging
from collections import OrderedDict
import json
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import cumtrapz
from scipy.interpolate import interp1d
from astropy.io import fits

from beagle_utils import prepare_data_saving, prepare_plot_saving, \
        BeagleDirectories, is_FITS_file, data_exists, plot_exists, set_plot_ticks


class BeagleMockCatalogue(object):

    def __init__(self, params_file):

        # Names of parameters, used to label the axes, whether they are log or
        # not, and possibly the extension name and column name containing the
        # "mock" parameters, i.e. the "true" parameters of the galaxy
        with open(params_file) as f:    
            # The use of the "OrderedDict" ensures that the order of the
            # entries in the dictionary reflects the order in the file
            self.adjust_params = json.load(f, object_pairs_hook=OrderedDict)

    def load(self, file_name=None):
        """ 
        Load a 'BEAGLE mock catalogue'

        Parameters
        ----------
        file_name : str
            Name of the file containing the catalogue.
        """

        if file_name is None:
            file_name = "BEAGLE_mock_catalogue.fits"

        logging.info("Loading the `BeagleMockCatalogue` file: " + file_name)

        name = file_name
        if not os.path.dirname(file_name):
            name = os.path.join(BeagleDirectories.results_dir, 
                    BeagleDirectories.pypbeagle_data, 
                    file_name)

        if is_FITS_file(name):
            self.data = fits.open(name)[1].data
            self.columns = fits.open(name)[1].columns
        else:
            self.data = ascii.read(name, Reader=ascii.basic.CommentedHeader)

    def compute(self, file_list, file_name=None, overwrite=False):
        """ 
        """ 

        if file_name is None:
            file_name = "BEAGLE_mock_catalogue.fits"

        # Check if the `file_name` already exists
        if data_exists(file_name) and not overwrite:
            logging.warning('The file `' + file_name + '` already exists. \n Exiting the function.')
            return

        # We save in a dictionary the extension name and column name containing
        # the "true" parameters
        params_dict = dict()
        # This cycles over all keys containing the parameter names
        for key, value in self.adjust_params.iteritems():
            d = { "extName" : "POSTERIOR PDF", "colName" : key}
            log = {"log" : False}
            log.update(value)
            # This cycles over the dictioary items for one parameter
            for in_key, in_value in value.iteritems():
                if in_key == "mock":
                    # This will merge the default dictionary `d` with the
                    # one found in the json file
                    d.update(in_value)

            params_dict[key] = d

        # Read all FITS file in the `file_list`, and extract the columns
        # defined in the `params_dict` dictionary
        data = OrderedDict()
        n_files = len(file_list)
        data['ID'] = np.chararray(n_files, itemsize=20)
        for i, file in enumerate(file_list):
            hdulist = fits.open(file)
            data['ID'][i] = os.path.basename(file).split('_BEAGLE')[0]
            for key, value in params_dict.iteritems():
                val = hdulist[value["extName"]].data[value["colName"]]
                if not key in data:
                    data[key] = np.zeros(n_files)
                data[key][i] = val
            hdulist.close()

        # Initialize a new (empty) primary HDU for your output FITS file
        hdulist = fits.HDUList(fits.PrimaryHDU())
    
        new_columns = list()
        for key, value in data.iteritems():

            # The `ID` column contains a string, while all the other columns real data
            if 'ID' in key:
                new_columns.append(fits.Column(name=str(key), format='20A', array=data[key]))
            else:
                new_columns.append(fits.Column(name=str(key), format='E', array=data[key]))

        cols_ = fits.ColDefs(new_columns)
        new_hdu = fits.BinTableHDU.from_columns(cols_)

        # And finally append the newly created binary table to the hdulist
        # that will be printed to the ouput FITS file
        hdulist.append(new_hdu)

        name = prepare_data_saving(file_name, overwrite=overwrite)
        hdulist.writeto(name, clobber=overwrite)

        self.columns = new_hdu.columns
        self.data = new_hdu.data

    def compare(self, summary_catalogue, 
            plot_name=None,
            summary_type='median',
            level=68.,
            params_to_plot=None,
            overwrite=False):

        # The summary statistics can be only 'mean' or 'median'
        if summary_type not in ('mean', 'median'):
            raise ValueError("`summary_type` con only be set to `mean` or `median`")

        if plot_name is None:
            plot_name = "BEAGLE_mock_retrieved_params.pdf"

        # Check if the plot already exists
        if plot_exists(plot_name) and not overwrite:
            logging.warning('The plot `' + plot_name + '` already exists. \n Exiting the function.')
            return

        # By default you plot all parameters
        if params_to_plot is None:
            params_to_plot = list()
            for key, value in self.adjust_params.iteritems():
                params_to_plot.append(key)

        _n = int(np.ceil(np.sqrt(len(params_to_plot))))
        fig, axs = plt.subplots(_n, _n)
        fig.subplots_adjust(left=0.08, bottom=0.08)
        fontsize = 8

        axs = np.ravel(axs)

        for i, param in enumerate(params_to_plot):
            # Extract the row corresponding to `param` from the mock catalogue
            # (this array contains the "true") values of the parameters

            true_param = self.data[param]
            _col = param + '_' + summary_type
            retrieved_param = summary_catalogue.hdulist['POSTERIOR PDF'].data[_col]

            _col = param + '_' + '{:.2f}'.format(level)
            error = summary_catalogue.hdulist['POSTERIOR PDF'].data[_col]

            ax = axs[i]

            # Set the x- and y-axis labels
            label = self.adjust_params[param]['label']
            xlabel = label + " (true)"
            ylabel = label

            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)

            # Check if we need to compute the log10 of the param or not
            x = true_param
            y = retrieved_param
            yerr_u = error[:,1]-y
            yerr_d = y-error[:,0]

            if 'log' in self.adjust_params[param]:
                if self.adjust_params[param]['log']:
                    x = np.log10(x)
                    y = np.log10(y)
                    yerr_d = np.log10(yerr_d)
                    yerr_u = np.log10(yerr_u)

            ax.errorbar(x, 
                    y, 
                    yerr=[yerr_d, yerr_u],
                    ls="",
                    marker='o',
                    ms=3, 
                    mew=0,
                    elinewidth=1.0,
                    capsize=3
                    )

            for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                         ax.get_xticklabels() + ax.get_yticklabels()):
                item.set_fontsize(fontsize)

            set_plot_ticks(ax, n_x=4, n_y=4)

        # Make the unused axes invisibles
        for ax in axs[len(params_to_plot):]:
            ax.axis('off')
        
        name = prepare_plot_saving(plot_name, overwrite=overwrite)

        plt.tight_layout()

        fig.savefig(name, dpi=None, facecolor='w', edgecolor='w',
                orientation='portrait', papertype='a4', format="pdf",
                transparent=False, bbox_inches="tight", pad_inches=0.1)

        plt.close(fig)
