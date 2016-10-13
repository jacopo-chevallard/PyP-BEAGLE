import os
import logging
from collections import OrderedDict
import json
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from astropy.io import fits
from  itertools import izip_longest

import bokeh.plotting as bk_plt
import bokeh.models as bk_mdl

from beagle_utils import prepare_data_saving, prepare_plot_saving, \
        BeagleDirectories, is_FITS_file, data_exists, plot_exists, set_plot_ticks

def grouper(n, iterable, fillvalue=None):
    "grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return izip_longest(fillvalue=fillvalue, *args)

def test_dim(testlist, dim=0):
   """tests if testlist is a list and how many dimensions it has
   returns -1 if it is no list at all, 0 if list is empty 
   and otherwise the dimensions of it"""
   if isinstance(testlist, list):
      if testlist == []:
          return dim
      dim = dim + 1
      dim = test_dim(testlist[0], dim)
      return dim
   else:
      if dim == 0:
          return -1
      else:
          return dim

def extract_err(err, data):
    '''private function to compute error bars
    Parameters
    ----------
    err : iterable
        xerr or yerr from errorbar
    data : iterable
        x or y from errorbar
    '''
    if (iterable(err) and len(err) == 2):
        a, b = err
        if iterable(a) and iterable(b):
            # using list comps rather than arrays to preserve units
            low = [thisx - thiserr for (thisx, thiserr)
                   in cbook.safezip(data, a)]
            high = [thisx + thiserr for (thisx, thiserr)
                    in cbook.safezip(data, b)]
            return low, high
    # Check if xerr is scalar or symmetric. Asymmetric is handled
    # above. This prevents Nx2 arrays from accidentally
    # being accepted, when the user meant the 2xN transpose.
    # special case for empty lists
    if len(err) > 1:
        fe = safe_first_element(err)
        if not ((len(err) == len(data) and not (iterable(fe) and
                                                len(fe) > 1))):
            raise ValueError("err must be a scalar, the same "
                             "dimensions as x, or 2xN.")
    # using list comps rather than arrays to preserve units
    low = [thisx - thiserr for (thisx, thiserr)
           in cbook.safezip(data, err)]
    high = [thisx + thiserr for (thisx, thiserr)
            in cbook.safezip(data, err)]
    return low, high

def errorbar(fig, x, y, xerr=None, yerr=None, kwargs={}):

  if xerr is not None:
      left, right = extract_err(xerr, x)
      fig.multi_line(left, right, **kwargs)

  if yerr is not None:
      lower, upper = extract_err(yerr, y)
      fig.multi_line(lower, upper, **kwargs)

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
            hdulist = fits.open(os.path.join(BeagleDirectories.results_dir, file))
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

    def compare_hist(self, summary_catalogue,
            bins=20,
            class_indices=None,
            class_colors=None,
            plot_name=None,
            summary_type='median',
            percentile=95.,
            params_to_plot=None,
            pow=False,
            overwrite=False,
            kwargs=None):

        # Check whether the IDs of the two catalogues match
        for i, ID in enumerate(self.data['ID']):
            if not ID == summary_catalogue.hdulist['POSTERIOR PDF'].data['ID'][i]:
                raise Exception("The object IDs of the `mock` and `summary` catalogues do not match!")

        # The summary statistics can be only 'mean' or 'median'
        if summary_type not in ('mean', 'median'):
            raise ValueError("`summary_type` con only be set to `mean` or `median`")

        if plot_name is None:
            plot_name = "BEAGLE_mock_retrieved_params_hist.pdf"

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
        fig.subplots_adjust(left=0.08, bottom=0.08, hspace=0.4)
        fontsize = 8

        axs = np.ravel(axs)

        diff_median = dict()
        diff_dispersion = dict()
        hist = list()
        edges = list()
        center = list()
        width = list()

        if kwargs is None:
            kwargs = {'alpha' : 0.7}

        for i, param in enumerate(params_to_plot):
            # Extract the row corresponding to `param` from the mock catalogue
            # (this array contains the "true") values of the parameters

            true_param = self.data[param]
            _col = param + '_' + summary_type
            retrieved_param = summary_catalogue.hdulist['POSTERIOR PDF'].data[_col]

            ax = axs[i]

            # Set the x- and y-axis labels
            label = self.adjust_params[param]['label']

            # Check if we need to compute the log10 of the param or not
            x = true_param
            y = retrieved_param

            if 'log' in self.adjust_params[param]:
                if self.adjust_params[param]['log'] and pow and '\log' in label:
                    label = label.replace('\log', '')
                if self.adjust_params[param]['log'] and not pow:
                    x = np.log10(x)
                    y = np.log10(y)

            xlabel = label + " (true-retrieved)"
            ylabel = '\# of objects'
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)

            diff = (x-y)

            if class_indices is None:
                class_indices = np.arange(len(diff))
                class_indices = [class_indices, ]

            diff_med = list()
            diff_disp = list()

            for i, indices in enumerate(class_indices):
                range = np.percentile(diff[indices], (0.5*(100.-percentile), 100.-0.5*(100.-percentile)))
                if class_colors is None:
                    color = None
                else:
                    color = class_colors[i]

                med, d0, d1 = np.percentile(diff[indices], [50., 0.5*(100.-68.), 100.-0.5*(100.-68.)])
                dispersion = 0.5*(d1-d0)

                diff_med.append(med)
                diff_disp.append(dispersion)

                ax.hist(diff[indices],
                    bins=bins,
                    range=range,
                    color=color,
                    lw='',
                    **kwargs)

                ax.text(0.5, 1.05, "$\sigma=" + "{:.3f}".format(dispersion) + "$", 
                    fontsize=8, 
                    horizontalalignment='center', 
                    transform=ax.transAxes)

                ax.text(0.5, 1.14, "$\\textnormal{med}=" + "{:.3f}".format(med) + "$", 
                    fontsize=8, 
                    horizontalalignment='center', 
                    transform=ax.transAxes)

            diff_median[param] = diff_med
            diff_dispersion[param] = diff_disp

            for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                         ax.get_xticklabels() + ax.get_yticklabels()):
                item.set_fontsize(fontsize)

            set_plot_ticks(ax, n_x=4, n_y=4)

        # Make the unused axes invisibles
        for ax in axs[len(params_to_plot):]:
            ax.axis('off')
        
        name = prepare_plot_saving(plot_name, overwrite=overwrite)

        #plt.tight_layout()

        fig.savefig(name, dpi=None, facecolor='w', edgecolor='w',
                orientation='portrait', papertype='a4', format="pdf",
                transparent=False, bbox_inches="tight", pad_inches=0.1)

        plt.close(fig)

        return diff_median, diff_dispersion

    def compare(self, summary_catalogue, 
            plot_name=None,
            summary_type='median',
            level=68.,
            params_to_plot=None,
            overwrite=False,
            rows=None, 
            interactive=False):

        # Check whether the IDs of the two catalogues match
        for i, ID in enumerate(self.data['ID']):
            if not ID == summary_catalogue.hdulist['POSTERIOR PDF'].data['ID'][i]:
                raise Exception("The object IDs of the `mock` and `summary` catalogues do not match!")

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

        # Do you consider only some rows in the catalogue?
        if rows is None:
            rows = np.arange(len(self.data.field(0)))

        _n = int(np.ceil(np.sqrt(len(params_to_plot))))

        ######################################################################
        # If interactive is `True`, you draw an interactive plot instead of
        # plotting into a pdf file
        ######################################################################
        if interactive:

            # Size (in pixels) of each Bokeh figure
            size = 400

            # Name of the output html file created by Bokeh
            name = prepare_plot_saving(os.path.splitext(plot_name)[0]+'.html', 
                    overwrite=overwrite)

            # Tell Bokeh to save the plot into an output html file 
            bk_plt.output_file(name)

            # create a column data source for the plots to share, and fill the
            # dictionary with all the data that will then be used in the Bokeh
            # plot
            data = dict()
            data['ID'] = self.data['ID'][rows]
            for param in params_to_plot:

                # Store the "true" parameter
                key = param + '_true'
                data[key] = self.data[param][rows]

                # and the retrieved one
                key = param + '_retrieved'
                _col = param + '_' + summary_type
                data[key] = summary_catalogue.hdulist['POSTERIOR PDF'].data[_col][rows]

                # Get the errors
                _col = param + '_' + '{:.2f}'.format(level)
                tmp = summary_catalogue.hdulist['POSTERIOR PDF'].data[_col][rows]

                # Store the errors in a way that is easily usable by the Bokeh
                # `multi_line` glyph
                err_xs = list()
                err_ys = list()
                key = param + '_true'
                for i, _x in enumerate(data[key]):
                    err_xs.append((_x, _x))
                    err_ys.append((tmp[i,0], tmp[i,1]))

                key = param + '_err_xs'
                data[key] = err_xs
                key = param + '_err_ys'
                data[key] = err_ys
                #data[key] = 0.5*(tmp[:,1]-tmp[:,0])

            strID = [ str(ID) + 
                    '_BEAGLE_input_for_Pierre.fits_snr_PS_CLEAR_PRISM_BEAGLE_marginal_SED_spec.pdf' 
                    for ID in data['ID']]

            data['strID'] = strID

            # Build the `ColumnDataSource` (a sort of dictionary) that will be
            # used by Bokeh to extract the data
            source = bk_mdl.ColumnDataSource(data=data)

            # Define the "tools" that will be included in the interactive Bokeh plot
            tools = 'wheel_zoom,pan,reset,resize,hover,tap'

            figs = list()
            for param in params_to_plot:
                # create a new plot and add a renderer
                fig = bk_plt.figure(tools=tools, width=size, height=size, title=None)

                # Plot the x and y data points as circles
                x = param + '_true'
                y = param + '_retrieved'
                error_low = param + '_err_xs'
                error_up = param + '_err_ys'
                fig.circle(x, y, source=source)

                # Plot the errorbars
                fig.multi_line(xs=error_low, ys=error_up, source=source, alpha=0.6)

                # Plot the 1-to-1 relation
                _max = np.amax(np.ravel([data[x], data[y]]))
                _min = np.amin(np.ravel([data[x], data[y]]))
                fig.line([_min, _max], [_min, _max], color='red')

                # Label the x- and y-axis
                label = self.adjust_params[param]['label_plain']

                xlabel = label + " (true)"
                fig.xaxis.axis_label = xlabel

                ylabel = label
                fig.yaxis.axis_label = ylabel
                    
                # Append the newly created figure to the `fgs` list of Bokeh figures
                figs.append(fig)

            # Arrange the different figures in a matrix-list
            grid_figs = list()

            for i in range(0,_n):
                _tmp = list()
                for j in range(_n):
                    if i*_n+j >= len(figs):
                        break
                    else:
                        _tmp.append(figs[i*_n+j])
                grid_figs.append(_tmp)

            # Use the matrix-list to create a grid of Bokeh figures
            p = bk_plt.gridplot(grid_figs)

            # Here we customize the behaviour of some tools
            hover = p.select(dict(type=bk_mdl.HoverTool))
            
            hover.tooltips = [
                 ("(x,y)", "($x, $y)"),
                ('ID', '@ID'),
            ]

            path = os.path.join(os.getcwd(),
                    '@strID')
            url = "file://" + path
            print "url: ", url
            taptool = p.select(type=bk_mdl.TapTool)
            taptool.callback = bk_mdl.OpenURL(url=url)

            bk_plt.show(p)

            return

        ######################################################################
        # Below is a standard plot written to a pdf file
        ######################################################################

        fig, axs = plt.subplots(_n, _n)
        fig.subplots_adjust(left=0.08, bottom=0.08)
        fontsize = 8

        axs = np.ravel(axs)

        for i, param in enumerate(params_to_plot):
            # Extract the row corresponding to `param` from the mock catalogue
            # (this array contains the "true") values of the parameters

            true_param = self.data[param][rows]
            _col = param + '_' + summary_type
            retrieved_param = summary_catalogue.hdulist['POSTERIOR PDF'].data[_col][rows]

            _col = param + '_' + '{:.2f}'.format(level)
            error = summary_catalogue.hdulist['POSTERIOR PDF'].data[_col][rows]

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
                    elinewidth=0.8,
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
