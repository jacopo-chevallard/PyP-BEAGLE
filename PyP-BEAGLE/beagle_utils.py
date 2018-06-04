import fnmatch
import re
import os
import logging
import numpy as np
from bisect import bisect_left
from scipy.integrate import simps, cumtrapz
from scipy.interpolate import interp1d
from datetime import datetime

import sys
import dependencies.WeightedKDE as WeightedKDE

import matplotlib.ticker as plticker


ID_COLUMN_LENGTH = 100

def is_integer(s):
    try:
        int(s)
        return True
    except ValueError:
        return False

def is_FITS_file(file_name):

    fits_ext = ('.fits', '.fits')
    arch_ext = ('', '.gz', '.z', '.zip')
    name = file_name.lower()
    for fExt in fits_ext:
        for aExt in arch_ext:
            suffix = fExt + aExt
            if name.endswith(suffix):
                return True

    return False

def trimFitsSuffix(file_name):

    fits_ext = ['.fits', '.fits']
    arch_ext = ['', '.gz', '.z', '.zip']
    name = file_name
    for fExt in fits_ext:
        for aExt in arch_ext:
            suffix = fExt + aExt
            if name.endswith(suffix):
                return name.split(suffix)[0]

    return False
            
def extract_IDs(table, key=None):

    ID_keys = ["ID", "ID_", "_ID"]

    if key is not None:
        return table[key]

    for name in table.dtype.names:
        if name == "ID" or name.startswith("ID") or name.endswith("ID"):
            return table[name]

    raise ValueError("Cannot extract a column containing the objects IDs!")

def extract_row(table, ID, key=None):

    if key is None:
        key='ID'

    if isinstance(table[key][0], basestring):
        row = table[table[key] == str(ID)]
    else:
        row = table[table[key] == int(ID)]

    if row is None or not row:
        raise ValueError("Cannot extract the row corresponding to the keyvalue '"+ID+"'")

    return row

def pause():

    try:
        wait = input("PRESS ENTER TO CONTINUE.")
    except KeyboardInterrupt:
        raise SystemExit

class BeagleDirectories(object):

    beagle_input_files = "BEAGLE-input-files"
    pypbeagle_data = os.path.join("pyp-beagle", "data")
    pypbeagle_plot = os.path.join("pyp-beagle", "plot")

    results_dir = ''

    suffix = 'BEAGLE'

    MN_suffix = '_BEAGLE_MNstats.dat'

    param_file = ''

    fontsize = 16
    inset_fontsize_fraction = 0.7

def get_files_list(results_dir=None, suffix=None):
    """ 
    Get all files ending with suffix.

    Parameters
    ----------
    results_dir : str, optional
        Directory containing the files to list. By default uses the
        RESULTS_DIR constant.

    suffix: str, optional
       Suffix of the files to list. Bu default ``BeagleDirectories.suffix``

    Returns
    -------

       file_list: str array
        List of all BEAGLE results files contained in the results directory.

       file_IDs: str array
        List of all IDs of corresponding to the BEAGLE results files contained
        in the results directory.
    """ 

    if results_dir is None:
        results_dir = BeagleDirectories.results_dir

    if suffix is None:
        suffix = BeagleDirectories.suffix + '.fits.gz'

    file_list = list()
    file_IDs = list()


    for file in sorted(os.listdir(results_dir)):
        if file.endswith(suffix) and os.path.getsize(os.path.join(results_dir, file)) > 0:
            file_list.append(file)
            file = file[0:file.find(suffix)-1]
            file_IDs.append(file)

    return file_list, file_IDs


def touch(file_name, times=None):
    with open(file_name, 'a'):
        os.utime(file_name, times)


def find_file(file_name, path):
    """ 
    Find a file in a directory.

    Parameters
    ----------
    file_name : str
        Name of the file to be found (can contain wildcards)

    path: str
        Path where to (recursively) search for file_name

    Returns
    -------
       : str
        Full path to the file_name
    """ 

    print "file_name: ", file_name

    for root, dirs, files in os.walk(path):
        for file in files:
            if fnmatch.fnmatch(file, file_name):
                return os.path.join(root, file)

    return None

def prepare_data_saving(file_name, results_dir=None, overwrite=False):
    """ 
    Prepare directory to save a data file.  

    Parameters
    ----------
    file_name : str
        Name of the output file (without directory tree).

    results_dir : str, optional
        Directory containing the BEAGLE output files. By default uses the
        RESULTS_DIR constant.

    overwrite: bool, optional
        If rue overwrites the file, is already present, while False makes a copy of
        the original file to avoid overwriting

    Returns
    -------
    name : str
        Full path to the output file,
    """ 

    if results_dir is None:
        results_dir = BeagleDirectories.results_dir

    directory = os.path.join(results_dir, BeagleDirectories.pypbeagle_data)
    if not os.path.exists(directory):
        logging.info("Creating the directory: " + directory)
        os.makedirs(directory)

    name = os.path.join(directory, os.path.basename(file_name))
    if os.path.isfile(name) and not overwrite:
        new_name = os.path.splitext(name)[0] + datetime.now().strftime("-%Y%m%d-%H%M%S") + os.path.splitext(name)[1]
        logging.warning("The file " + name + " already exists, and it will be renamed to " + new_name)
        os.rename(name, new_name)

    return name

def getPathForPlot(file_name, results_dir=None):
    """ 
    """ 

    if results_dir is None:
        results_dir = BeagleDirectories.results_dir

    directory = os.path.join(results_dir, BeagleDirectories.pypbeagle_plot)
    return os.path.join(directory, os.path.basename(file_name))

def plot_exists(file_name, results_dir=None):
    """ 
    Check if a plot already exists in the PyP-BEAGLE tree.

    Parameters
    ----------
    file_name : str
        Name of the output file (without directory tree).

    results_dir : str, optional
        Directory containing the BEAGLE output files. By default uses the
        RESULTS_DIR constant.

    Returns
    -------
    bool
        Whether the plot aready exists or not.
    """ 

    return os.path.isfile(getPathForPlot(file_name, results_dir))

def getPathForData(file_name, results_dir=None):
    """ 
    """ 

    if results_dir is None:
        results_dir = BeagleDirectories.results_dir

    directory = os.path.join(results_dir, BeagleDirectories.pypbeagle_data)
    return os.path.join(directory, os.path.basename(file_name))

def data_exists(file_name, results_dir=None):
    """ 
    Check if a file already exists in the PyP-BEAGLE tree.

    Parameters
    ----------
    file_name : str
        Name of the output file (without directory tree).

    results_dir : str, optional
        Directory containing the BEAGLE output files. By default uses the
        RESULTS_DIR constant.

    Returns
    -------
    bool
        Whether the plot aready exists or not.
    """ 

    return os.path.isfile(getPathForData(file_name, results_dir))

def prepare_plot_saving(file_name, results_dir=None, overwrite=False):
    """ 
    Prepare directory to save a plot file.  

    Parameters
    ----------
    file_name : str
        Name of the output file (without directory tree).

    results_dir : str, optional
        Directory containing the BEAGLE output files. By default uses the
        RESULTS_DIR constant.

    overwrite: bool, optional
        If rue overwrites the file, is already present, while False makes a copy of
        the original file to avoid overwriting

    Returns
    -------
    name : str
        Full path to the output file,
    """ 

    if results_dir is None:
        results_dir = BeagleDirectories.results_dir

    directory = os.path.join(results_dir, BeagleDirectories.pypbeagle_plot)
    if not os.path.exists(directory):
        logging.info("Creating the directory: " + directory)
        os.makedirs(directory)

    name = os.path.join(directory, os.path.basename(file_name))
    if os.path.isfile(name) and not overwrite:
        new_name = os.path.splitext(name)[0] + datetime.now().strftime("-%Y%m%d-%H%M%S") + os.path.splitext(name)[1]
        logging.warning("The file " + name + " already exists, and it will be renamed to " + new_name)
        os.rename(name, new_name)

    return name

def set_font_size(ax, fontsize):

    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                 ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(fontsize)


def set_plot_ticks(ax, which='both', n_x=4, n_y=4, prune_x=None, prune_y=None):
    """ 
    Set the major and minor ticks of plots.

    Parameters
    ---------
    ax : Axes instance
        The axis on which to apply yhe major/minor tick marks setting.

    n_x : int, optional
        Maximum number of intervals on the x-axis.

    n_y : int, optional
        Maximum number of intervals on the y-axis.

    prune_x : str, optional 
        Either 'lower' | 'upper' | 'both' | None, to remove
        edge ticks on the x-axis.

    prune_y : str, optional 
        Either 'lower' | 'upper' | 'both' | None, to remove
        edge ticks on the y-axis.
    """ 

    if which.lower() == 'both' or which.lower() == 'x':
        ax.xaxis.major_locations = plticker.MaxNLocator(nbins=n_x, prune=prune_x) 
        ax.xaxis.set_major_locator(ax.xaxis.major_locations)

        ax.xaxis.minor_locations = plticker.AutoMinorLocator(2)
        ax.xaxis.set_minor_locator(ax.xaxis.minor_locations)

    if which.lower() == 'both' or which.lower() == 'y':
        ax.yaxis.major_locations = plticker.MaxNLocator(nbins=n_y, prune=prune_y) 
        ax.yaxis.set_major_locator(ax.yaxis.major_locations)

        ax.yaxis.minor_locations = plticker.AutoMinorLocator(2)
        ax.yaxis.set_minor_locator(ax.yaxis.minor_locations)

def match_ID(ID_list_1, ID_list_2, sorted=False, ignore_string=None):
    """ 
    Match the ID in two catalogues.

    Parameters
    ----------
    ID_list_1 : numpy array int
        The ID for the first catalogue.

    ID_list_2 : numpy array int
        The ID for the second catalogue.

    sorted : bool, or list of bool, optional
        If True then the ID arrays are assumed to be already sorted, otherwise
        they are sorted in the routine 

    Returns
    -------
    indices_1 : numpy array int
        Array of indices such as `ID_list_1`[indices_1] = `ID_list_2`[indices_2]

    indices_2 : numpy array int
        Array of indices such as `ID_list_1`[indices_1] = `ID_list_2`[indices_2]

    """

    # Firstly, check weather ID_list_1 is longer than ID_list_2 or viceversa
    n1 = len(ID_list_1)
    n2 = len(ID_list_2)

    if n1 >= n2:
        ID_long = np.array(ID_list_1)
        n_long = n1
        ID_short = np.array(ID_list_2)
        n_short = n2
    else:
        ID_long = np.array(ID_list_2)
        n_long = n2
        ID_short = np.array(ID_list_1)
        n_short = n1

    is_int = False
    sort_long = range(n_long)
    sort_short = range(n_short)
    if issubclass(ID_long.dtype.type, np.integer) and issubclass(ID_short.dtype.type, np.integer):
        is_int = True
        if not sorted:
            # Sort the two arrays and get the indices
            sort_long = np.argsort(ID_long)
            sort_short = np.argsort(ID_short)
    else:
        ID_long = np.array(ID_long, dtype=str)
        ID_short = np.array(ID_short, dtype=str)

    match_indx_long = np.full(n_long, -1, dtype=np.int)
    match_indx_short = np.full(n_short, -1, dtype=np.int)
    mask_long = np.zeros(n_long, dtype=bool)

    for i in range(n_short):

        if is_int:
            i1 = bisect_left(ID_long[sort_long], ID_short[sort_short[i]])

            if i1 < n_long:
                if ID_long[sort_long[i1]] == ID_short[sort_short[i]]:
                    match_indx_long[i1] = sort_long[i1]
                    match_indx_short[i] = sort_short[i]
        else:
            for j in range(n_long):
                if mask_long[j]:
                    continue

                if ignore_string is not None:
                    if isinstance(ignore_string, re.compile('').__class__):
                        ID_l = ignore_string.sub('', ID_long[j])
                        ID_s = ignore_string.sub('', ID_short[i])
                    else:
                        ID_l = ID_long[j].replace(ignore_string, '')
                        ID_s = ID_short[i].replace(ignore_string, '')
                else:
                    ID_l = ID_long[j]
                    ID_s = ID_short[i]

                if ID_l == ID_s:
                    match_indx_long[i] = j
                    match_indx_short[i] = i
                    mask_long[j] = True

    if n1 >= n2:
        indices_1 = match_indx_long[match_indx_long >= 0]
        indices_2 = match_indx_short[match_indx_short >= 0]
    else:
        indices_1 = match_indx_short[match_indx_short >= 0]
        indices_2 = match_indx_long[match_indx_long >= 0]

    return indices_1, indices_2


def weighted_avg_and_std(values, weights):
    """
    Return the weighted average and standard deviation.

    Parameters
    ----------
    values : Numpy ndarray
        Contains the value of the parameter.

    weights : Numpy ndarray
        Contains the weights. Must have same shape as values.

    Returns
    -------
    average : float 
        The weighted average of the input values.

    stddev : float
        The weighted standard deviation of the input values.
    """

    average = np.average(values, weights=weights)
    variance = np.average((values-average)**2, weights=weights)  # Fast and numerically precise

    return (average, np.sqrt(variance))


def prepare_violin_plot(data, 
        weights=None, 
        min_x=None,
        max_x=None,
        nXgrid=100,
        max_interval=99.7):

    if min_x is None:
        min_x = np.min(data) 

    if max_x is None:
        max_x = np.max(data) 

    # Compute the marginal PDF through a weighted KDE
    if weights is not None:
        pdf = WeightedKDE.gaussian_kde(data, weights=weights)
    else:
        pdf = WeightedKDE.gaussian_kde(data)

    # Build a grid of value over which computing the actual PDF from its KDE
    x_grid = np.linspace(min_x, max_x, nXgrid)

    # Compute the PDF
    pdf_grid = np.array(pdf(x_grid))

    # ******************************************************************
    # NB: in this case it is correct to integrate the pdf, and not just to sum,
    # since the result of the weighted KDE algorithm is a *density*. You can
    # check this by checking the value of `pdf_norm` below, which is indeed ~1.
    # This also justify the use of `cumtrapz` below (and not just of `cumsum`,
    # as in `beagle_summary_catalogue.get1DInterval`).
    # ******************************************************************

    # Compute the PDF normalization, so its integral will be exactly 1
    pdf_norm = simps(pdf_grid, x_grid)
    pdf_grid /= pdf_norm

    # Now compute the cumulative PDF
    cumul_pdf = cumtrapz(pdf_grid, x_grid, initial=0.)
    cumul_pdf /= cumul_pdf[-1]

    # Get the interpolant of the cumulative PDF
    f_interp = interp1d(cumul_pdf, x_grid)

    # Compute the limits over which you will plot the "violin", for
    # instance showing only the cumulative PDF up to +/- 3 sigma
    intv = 0.5*(100.-max_interval)/100.
    lims = f_interp([intv,1.-intv])
    prob_lims = pdf(lims) / pdf_norm

    # The median corresponds to a cumulative probability = 0.5
    median = f_interp(0.5)

    i1 = bisect_left(x_grid, lims[0])
    i2 = bisect_left(x_grid, lims[1])

    x_violin = np.concatenate(([lims[0]], x_grid[i1+1:i2], [lims[1]]))

    y_violin = np.concatenate(([prob_lims[0]], pdf_grid[i1+1:i2], [prob_lims[1]]))


    return pdf, pdf_norm, median, x_violin, y_violin
