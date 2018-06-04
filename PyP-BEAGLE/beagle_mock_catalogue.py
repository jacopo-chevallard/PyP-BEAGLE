import os
import logging
from collections import OrderedDict
import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from scipy.interpolate import interp1d
from astropy.io import fits
from  itertools import izip_longest

import bokeh.plotting as bk_plt
import bokeh.models as bk_mdl

from beagle_utils import prepare_data_saving, prepare_plot_saving, \
        BeagleDirectories, is_FITS_file, data_exists, plot_exists, set_plot_ticks, \
        is_integer, match_ID, ID_COLUMN_LENGTH

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

    def __init__(self, params_file, 
            ID_key="ID", 
            file_name=None,
            ignore_string=None, 
            overwrite_plots=True,
            plot_title=None,
            n_bins=10):

        # Names of parameters, used to label the axes, whether they are log or
        # not, and possibly the extension name and column name containing the
        # "mock" parameters, i.e. the "true" parameters of the galaxy
        with open(params_file) as f:    
            # The use of the "OrderedDict" ensures that the order of the
            # entries in the dictionary reflects the order in the file
            self.adjust_params = json.load(f, object_pairs_hook=OrderedDict)

        self.ID_key = ID_key

        self.ignore_string = ignore_string

        self.overwrite_plots = overwrite_plots

        self.n_bins = n_bins

        self.plot_title = None
        if plot_title is not None:
            self.plot_title = plot_title
            if rcParams['text.usetex']:
                self.plot_title = self.plot_title.replace('_', '\_')

        if file_name is not None:
            self.load(file_name)

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

        self.file_name = os.path.basename(file_name)

        logging.info("Loading the `BeagleMockCatalogue` file: " + file_name)

        name = file_name
        if not os.path.dirname(file_name):
            name = os.path.join(BeagleDirectories.results_dir, 
                    BeagleDirectories.pypbeagle_data, 
                    file_name)

        self.hdulist = None

        if is_FITS_file(name):
            hdulist = fits.open(name)
            if len(hdulist) > 2:
                self.hdulist = hdulist
            else:
                self.data = fits.open(name)[1].data
                self.columns = fits.open(name)[1].columns
        else:
            self.data = ascii.read(name, Reader=ascii.basic.CommentedHeader)

    def get_param_values(self, ID, names):

        values = np.zeros(len(names), dtype=np.float32)

        for i, name in enumerate(names):

            if "extName" in self.adjust_params[name]:
                extName = self.adjust_params[name]["extName"]
            else:
                extName = "POSTERIOR PDF"

            if "colName" in self.adjust_params[name]:
                colName = self.adjust_params[name]["colName"]
            else:
                colName = name


            if self.hdulist is not None:
                data = self.hdulist[extName].data
            else:
                data = self.data

            if self.ID_key in data.dtype.names:
                if self.ignore_string is not None:
                    values[i] = data[colName][data[self.ID_key] == ID.replace(self.ignore_string, '')]
                else:
                    values[i] = data[colName][data[self.ID_key] == ID]
            else:
                if is_integer(ID):
                    row = int(ID)
                else:
                    row = int(ID.split('_')[0])

                values[i] = data[colName][row-1]

        return values

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
                new_columns.append(fits.Column(name=str(key), format=str(ID_COLUMN_LENGTH)+'A', array=data[key]))
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
            class_indices=None,
            class_colors=None,
            plot_name=None,
            summary_type='median',
            percentile=95.,
            params_to_plot=None,
            pow=False,
            overwrite=False,
            kwargs=None):
            

        # Match IDs in the two catalogs
        indx1, indx2  = match_ID(self.hdulist[1].data['ID'], summary_catalogue.hdulist[1].data['ID'], ignore_string=self.ignore_string)

        #print "indx1: ", indx1, len(self.hdulist[1].data['ID']), max(indx1), len(indx1)
        #print "indx2: ", indx2, len(summary_catalogue.hdulist[1].data['ID']), max(indx2), len(indx2)

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

            if "extName" in self.adjust_params[param]:
                extName = self.adjust_params[param]["extName"]
            else:
                extName = "POSTERIOR PDF"

            if "colName" in self.adjust_params[param]:
                colName = self.adjust_params[param]["colName"]
            else:
                colName = param

            true_param = self.hdulist[extName].data[colName][indx1]

            _col = colName + '_' + summary_type

            retrieved_param = summary_catalogue.hdulist[extName].data[_col][indx2]

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
                ind_set = set(indices)
                #print "indices: ", indices, len(diff)
                idx = [j for j, item in enumerate(indx1) if item in ind_set]
                indx = idx
                #print "indx: ", indx

                if 'diff_range' in self.adjust_params[param]:
                    range = self.adjust_params[param]['diff_range']
                else:
                    range = np.percentile(diff[idx], (0.5*(100.-percentile), 100.-0.5*(100.-percentile)))

                if class_colors is None:
                    color = None
                else:
                    color = class_colors[i]

                med, d0, d1 = np.percentile(diff[idx], [50., 0.5*(100.-68.), 100.-0.5*(100.-68.)])
                dispersion = 0.5*(d1-d0)

                diff_med.append(med)
                diff_disp.append(dispersion)

                ax.hist(diff[idx],
                    bins=self.n_bins,
                    range=range,
                    color=color,
                    lw=0,
                    **kwargs)

                ax.text(0.1+i*0.25, 1.05, "$\sigma=" + "{:.3f}".format(dispersion) + "$", 
                    fontsize=8, 
                    horizontalalignment='center', 
                    transform=ax.transAxes)

                ax.text(0.1+i*0.25, 1.14, "$\\textnormal{med}=" + "{:.3f}".format(med) + "$", 
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
            

        # Match IDs in the two catalogs
        indx1, indx2  = match_ID(self.hdulist[1].data['ID'], summary_catalogue.hdulist[1].data['ID'], ignore_string=self.ignore_string)

        # Check whether the IDs of the two catalogues match
        #if 'ID' in self.hdulist[1].data:
        #    for i, ID in enumerate(self.hdulist[1].data['ID']):
        #        if not ID == summary_catalogue.hdulist['POSTERIOR PDF'].data['ID'][i]:
        #            raise Exception("The object IDs of the `mock` and `summary` catalogues do not match!")

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
        #if rows is None:
        #    rows = np.arange(len(self.hdulist[1].data.field(0)))

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
            data['ID'] = self.hdulist[1].data['ID'][indx1]
            for param in params_to_plot:

                # Store the "true" parameter
                key = param + '_true'
                data[key] = self.data[param][indx1]

                # and the retrieved one
                key = param + '_retrieved'
                _col = param + '_' + summary_type
                data[key] = summary_catalogue.hdulist['POSTERIOR PDF'].data[_col][indx2]

                # Get the errors
                _col = param + '_' + '{:.2f}'.format(level)
                tmp = summary_catalogue.hdulist['POSTERIOR PDF'].data[_col][indx2]

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

            if "extName" in self.adjust_params[param]:
                extName = self.adjust_params[param]["extName"]
            else:
                extName = "POSTERIOR PDF"

            if "colName" in self.adjust_params[param]:
                colName = self.adjust_params[param]["colName"]
            else:
                colName = param

            true_param = self.hdulist[extName].data[colName][indx1]

            _col = colName + '_' + summary_type

            retrieved_param = summary_catalogue.hdulist[extName].data[_col][indx2]

            _col = colName + '_' + '{:.2f}'.format(level)
            error = summary_catalogue.hdulist[extName].data[_col][indx2]

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

            is_log = False
            if 'log' in self.adjust_params[param]:
                if self.adjust_params[param]['log']:
                    is_log = True
                    #x = np.log10(x)
                    #y = np.log10(y)
                    #yerr_d = np.log10(yerr_d)
                    #yerr_u = np.log10(yerr_u)

            ax.errorbar(x, 
                    y, 
                    yerr=[yerr_d, yerr_u],
                    ls="",
                    marker='o',
                    ms=3, 
                    mew=0,
                    elinewidth=0.6,
                    capsize=3,
                    alpha=0.5
                    )

            for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                         ax.get_xticklabels() + ax.get_yticklabels()):
                item.set_fontsize(fontsize)

            # Make the axes in log-scale
            if is_log:
                ax.set_xscale('log')
                ax.set_yscale('log')
            else:
                set_plot_ticks(ax, n_x=4, n_y=4)

            # Make axes of same aspect ratio
            plt.axis('equal')

            # Set same axis range for x- and y-axis, and use the (optional)
            # user-provided value  
            xx = ax.get_xlim()
            yy = ax.get_ylim()
            ll = [min([xx[0], yy[0]]), max(xx[1], yy[1])]
            if 'axis_range' in self.adjust_params[param]:
                rr = self.adjust_params[param]['axis_range']
                ll = [max(ll[0],rr[0]), min(ll[1],rr[1])]
                ax.set_autoscale_on(False)

            ax.set_xlim(ll)
            ax.set_ylim(ll)

            # Plot the 1-to-1 relations
            ax.plot(ax.get_xlim(), ax.get_ylim(), 
                    ls="-",
                    color="black",
                    c=".3")

        # Make the unused axes invisibles
        for ax in axs[len(params_to_plot):]:
            ax.axis('off')
        
        name = prepare_plot_saving(plot_name, overwrite=overwrite)

        plt.tight_layout()

        fig.savefig(name, dpi=None, facecolor='w', edgecolor='w',
                orientation='portrait', papertype='a4', format="pdf",
                transparent=False, bbox_inches="tight", pad_inches=0.1)

        plt.close(fig)


    def plot_input_param_distribution(self, savefig=True):

        n_par = len(self.adjust_params)
        fig, axs = plt.subplots(n_par, n_par)
        fig.subplots_adjust(left=0.08, bottom=0.08, hspace=0.2)
        fontsize = 8
        axes_linewidth = 0.7
        color = "grey"

        last = n_par-1

        for i, (keyX, valueX) in enumerate(self.adjust_params.iteritems()):

            if "extName" in self.adjust_params[keyX]:
                extNameX = self.adjust_params[keyX]["extName"]
            else:
                extNameX = "POSTERIOR PDF"

            if "colName" in self.adjust_params[keyX]:
                colNameX = self.adjust_params[keyX]["colName"]
            else:
                colNameX = keyX

            valuesX = self.hdulist[extNameX].data[colNameX]
            if 'log' in self.adjust_params[keyX]:
                if self.adjust_params[keyX]['log']:
                    valuesX = np.log10(valuesX)

            xlabel = self.adjust_params[keyX]['label']

            for j, (keyY, valueY) in enumerate(self.adjust_params.iteritems()):

                if "extName" in self.adjust_params[keyY]:
                    extNameY = self.adjust_params[keyY]["extName"]
                else:
                    extNameY = "POSTERIOR PDF"

                if "colName" in self.adjust_params[keyY]:
                    colNameY = self.adjust_params[keyY]["colName"]
                else:
                    colNameY = keyY

                valuesY = self.hdulist[extNameY].data[colNameY]
                if 'log' in self.adjust_params[keyY]:
                    if self.adjust_params[keyY]['log']:
                        valuesY = np.log10(valuesY)

                ylabel = self.adjust_params[keyY]['label']

                ax = axs[i][j]

                if j > i:
                    ax.axis('off')
                    continue

                elif i == j:

                    ax.hist(valuesX,
                        bins=self.n_bins,
                        lw=0,
                        color=color)

                elif j < i:

                    ax.plot(valuesY,
                            valuesX,
                            marker='o',
                            ms=2, 
                            color=color,
                            ls="")


                    #ax.set_xlim(rangeY)
                    #ax.set_ylim(rangeX)

                if i == j:
                    ax.yaxis.tick_right()
                    if i == 0:
                        ax.set_xticklabels(('',))
                        ax.set_ylabel(xlabel)
                    elif i == last:
                        ax.set_xlabel(ylabel)
                    else:
                        ax.set_xticklabels(('',))
                elif j == 0:
                    ax.set_ylabel(xlabel)
                    if i != last:
                        ax.set_xticklabels(('',))
                elif i == last:
                    ax.set_xlabel(ylabel)
                    ax.set_yticklabels(('',))
                else:
                    ax.set_xticklabels(('',))
                    ax.set_yticklabels(('',))


                #if i == last:
                #    ax.set_xlabel(ylabel)
                #    if j != 0 and j != last:
                #        ax.set_yticklabels(('',))
                #    elif j == last:
                #        ax.yaxis.tick_right()

                #ax.set_ylabel('j ' + str(j) + 'i ' + str(i) )
                #if i == 0:

                #if j != 0 and i != (n_par-1) and i != j:
                #    ax.set_xticklabels(('',))
                #    ax.set_yticklabels(('',))

                #if i == j and i != 0 and i != (n_par-1):
                #    if j != 0:
                #        ax.set_xticklabels(('',))

                #if j == 0:
                #    ax.set_ylabel(xlabel)
                #    if i != 0:
                #            ax.set_xticklabels(('',))


                for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                             ax.get_xticklabels() + ax.get_yticklabels()):
                    item.set_fontsize(fontsize)

                set_plot_ticks(ax, n_x=3, n_y=3)

                ax.tick_params(which='minor', axis='both',
                                length=1.5, width=axes_linewidth)

                ax.tick_params(which='major', axis='both',
                                length=3, width=axes_linewidth)

                for axis in ['top','bottom','left','right']:
                  ax.spines[axis].set_linewidth(axes_linewidth)


        # Make the unused axes invisibles
        #for ax in axs[n_par:]:
        #    ax.axis('off')

        # Fine-tune figure; hide x ticks for top plots and y ticks for right plots
        plt.setp([a.get_xticklabels() for a in axs[0, :]], visible=False)
        plt.setp([a.get_yticklabels() for a in axs[:, 1]], visible=False)
    
        if self.plot_title is not None:
            fig.suptitle(self.plot_title)

        if savefig:

            plot_name = "BEAGLE_mock_input_params_distribution.pdf"
            name = prepare_plot_saving(plot_name, overwrite=self.overwrite_plots)

            #plt.tight_layout()

            fig.savefig(name, dpi=None, facecolor='w', edgecolor='w',
                    orientation='portrait', papertype='a4', format="pdf",
                    transparent=False, bbox_inches="tight", pad_inches=0.1)

        plt.close(fig)

        return fig, axs
