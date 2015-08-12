import os
import logging
import numpy as np
from bisect import bisect_left
from datetime import datetime


def prepare_file_writing(results_dir, file_name):
    """ 
    Prepare directory for writing the file.  

    Parameters
    ----------
    results_dir : str
        Directory containing the BANGS output files.

    file_name : str
        Name of the output file (without directory tree).

    Returns
    -------
    name : str
        Full path to the output file,
    """ 

    directory = os.path.join(results_dir, 'pybangs', 'data')
    if not os.path.exists(directory):
        logging.info("Creating the directory: " + directory)
        os.makedirs(directory)

    name = os.path.join(directory, os.path.basename(file_name))
    if os.path.isfile(name):
        new_name = os.path.splitext(name)[0] + datetime.now().strftime("-%Y%m%d-%H%M%S") + os.path.splitext(name)[1]
        logging.warning("The file " + name + " already exists, and it will be renamed to " + new_name)
        os.rename(name, new_name)

    return name

def match_ID(ID_list_1, ID_list_2, sorted=False ):

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
        ID_long = ID_list_1
        n_long = n1
        ID_short = ID_list_2
        n_short = n2
    else:
        ID_long = ID_list_2
        n_long = n2
        ID_short = ID_list_1
        n_short = n1

    if not sorted:
        # Sort the two arrays and get the indices
        sort_long = np.argsort(ID_long)
        sort_short = np.argsort(ID_short)
    else:
        sort_long = range(n_long)
        sort_short = range(n_short)

    match_indx_long = np.full(n_long, -1, dtype=np.int)
    match_indx_short = np.full(n_short, -1, dtype=np.int)

    for i in range(n_short):
        i1 = bisect_left(ID_long[sort_long], ID_short[sort_short[i]])

        if ID_long[sort_long[i1]] == ID_short[sort_short[i]]:
            match_indx_long = sort_long[i1]
            match_indx_short = sort_short[i]

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
