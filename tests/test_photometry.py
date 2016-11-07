import sys
sys.path.append("../PyP-BEAGLE")
import os
import argparse
import glob
import logging
import ConfigParser
import multiprocessing as mp
import copy_reg
import types
from matplotlib import rc
from astropy.io import ascii
from astropy.io import fits
import numpy as np
import random

from beagle_photometry import Photometry
from beagle_pdf import PDF
from beagle_utils import BeagleDirectories, get_files_list
from beagle_parsers import standard_parser

def _reduce_method(m):
    if m.im_self is None:
        return getattr, (m.im_class, m.im_func.func_name)
    else:
        return getattr, (m.im_self, m.im_func.func_name)

if __name__ == '__main__':

    parser = standard_parser()

    # Get parsed arguments
    args = parser.parse_args()    

    # Set logging level
    logging.basicConfig(level=args.loglevel)

    # Set directory containing BEAGLE results files
    BeagleDirectories.results_dir = args.results_dir

    # Read parameter file
    config = ConfigParser.SafeConfigParser()
    
    # Save name of the parameter file
    param_file = args.param_file
    BeagleDirectories.param_file = param_file 

    # Create the full path to the parameter file, taking advantage of the fact
    # that BEAGLE makes a copy of the input parameter file in the results
    # folder
    name = os.path.join(BeagleDirectories.results_dir, BeagleDirectories.beagle_input_files, param_file) 

    # Parse the parameter file
    # NB: remember that to access the parsed parameter file, the first line of the file must contain the string
    # [main]
    # without the "#" symbol
    config.read(name)

    # Set global font size
    font = {'size': 16}
    rc('font', **font)

    # Get list of results files and object IDs from the results directory (it
    # automatically exclude files with size 0)
    file_list, IDs = get_files_list()

    # Initialize an instance of the main "Photometry" class
    my_photometry = Photometry()

    # Set parameter names and labels
    my_PDF = PDF(os.path.join(BeagleDirectories.results_dir, "params_names.json"))

    # Read name of filter file from the parsed parameter file
    filters_file = os.path.expandvars(config.get('main', 'FILTERS FILE'))

    # Load a set of photometric filters
    my_photometry.filters.load(filters_file)

    # Read name of observed catalogue from the parsed parameter file
    file_name = os.path.expandvars(config.get('main', 'PHOTOMETRIC CATALOGUE'))

    # Load observed catalogue
    my_photometry.observed_catalogue.load(file_name)

    # Cycle across all files found in the results directory and create the
    # triangle and marginal plots
    for ID in IDs:
        my_photometry.plot_marginal(ID)
        my_PDF.plot_triangle(ID, M_star=True)

    # Uncomment the following line if you need to perform further operations
    # with PyP-BEAGLE (see below)
    sys.exit()

    # Compute the "summary catalogue", i.e. summary statistics for all the
    # quantities deifned in the "summary_config.json" file
    my_photometry.summary_catalogue.compute(file_list)

    # In some occasions you may just want to load a pre-computed summary catalogue
    #my_photometry.summary_catalogue.load(file_name)

    #########################################
    ############### MultiNest catalogue - start
    #########################################
    # This part is specific for analyses for which one need to consider the
    # presence of multi-modal solutions
    # Contact the PyP-BEAGLE authors for details

    #file_name = "UVUDF_MultiNest.cat"

    #my_photometry.multinest_catalogue.load(file_name)

    # ********* Compute ***************
    #file_list = list()
    #for file in os.listdir(results_dir):
    #    if file.endswith("MNstats.dat"):
    #        file_list.append(file)

    #n_par = 6
    #my_photometry.multinest_catalogue.compute( n_par, file_list, file_name)

    #########################################
    ############### MultiNest catalogue - end
    #########################################

    # Compute posterior predictive checks (see Sec. 5.1 of arxiv.org/abs/1603.03037)
    my_photometry.PPC.compute(my_photometry.observed_catalogue, 
            my_photometry.filters)

    # Load pre-computed posterior predictive checks
    #my_photometry.PPC.load( file_name) 

    # Plot the chi-square test statistics
    my_photometry.PPC.plot_chi2()

    # Plot the integrated tail probability (p-value)
    my_photometry.PPC.plot_p_value(broken_axis=True)

    # Compute residual photometry, i.e. difference between data and simulated
    # data (see Sec. 5.1 of arxiv.org/abs/1603.03037)
    my_photometry.residual.compute(my_photometry.observed_catalogue,
            my_photometry.summary_catalogue, my_photometry.filters)
