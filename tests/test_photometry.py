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

    param_file = args.param_file
    BeagleDirectories.param_file = param_file 

    print "param_file: ", param_file
    name = os.path.join(BeagleDirectories.results_dir, param_file) 
    config.read(name)

    # Global font size
    font = {'size': 16}
    rc('font', **font)

    # Get list of results files and object IDs from the results directory
    file_list, IDs = get_files_list()

    # Initialize an instance of the main "Photometry" class
    my_photometry = Photometry()

    # Set parameter names and labels
    my_PDF = PDF(os.path.join(BeagleDirectories.results_dir, "params_names.json"))

    # We now have access to other classes!

    # *****************************************************
    # *********** Filters ****************
    # *****************************************************

    # We can load a set of photometric filters
    filters_file = os.path.expandvars(config.get('main', 'FILTERS FILE'))

    my_photometry.filters.load(filters_file)

    # *****************************************************
    # *********** Observed Catalogue ****************
    # *****************************************************

    # ********** Loading *****************
    file_name = os.path.expandvars(config.get('main', 'PHOTOMETRIC CATALOGUE'))
    print "file_name: ", file_name
    my_photometry.observed_catalogue.load(file_name)

   # random.seed(12345)
   # rand_IDs = random.sample(IDs, 200)
   # rand_IDs = ('XDFB-3345164827', 'XDFB-3701762974', 'XDFI-3751560108', 'XDFI-3827975128', 
   #         'XDFV-3676055820', 'XDFV-4264462281', 'XDFY-4313762847', 'XDFY-4307762415', 
   #         'XDFZ-4256373146', 'XDFZ-3870255364')

    rand_IDs = IDs
    #rand_IDs = ('COSMOS', )

    specFiles = glob.glob("/Users/jchevall/Coding/BEAGLE/files/results/JWST/March_2016/Pierre/*.fits.gz")

    for ID in rand_IDs:
        if any(ID in spec for spec in specFiles):
            my_photometry.plot_marginal(ID)
            my_PDF.plot_triangle(ID, M_star=True)
    stop

    #pool = mp.Pool(4)
    #results = list()

    #copy_reg.pickle(types.MethodType, _reduce_method)

    #for ID in IDs[0:7]:
       # print "ID: ", ID
        #pool.apply_async(print_ID, args=[config])
        #pool.map(getattr_proxy(my_photometry, "print_ID"), ID)
        #pool.apply_async(getattr_proxy(my_photometry, "plot_marginal"), ID)
        #pool.apply_async(my_photometry.plot_marginal, args=[my_photometry.plot_marginal, ID])
        #pool.apply_async(wrap_make_photometry_plots, (ID, my_photometry, my_PDF))

        #pool.apply_async(my_print_ID, (ID,), dict(photometry=my_photometry, PDF=my_PDF))
     #   my_other_print_ID(ID, my_photometry)
        #pool.apply_async(my_other_print_ID, args=[ID, my_photometry], callback=results.append)

    #while len(results) < 8:
    #    continue

    #pool.close()

    # ********** Plotting of the marginal photometry *****************
    #my_photometry.plot_marginal(ID)
    #my_photometry.plot_replicated_data(ID)

    # *****************************************************
    # *********** "BANGS summary catalogue" ****************
    # *****************************************************

    file_name = "BANGS_summary_catalogue.fits"

    # ********* Load ***************
    #my_photometry.summary_catalogue.load(file_name)

    # ********* Compute ***************
    #file_list = ("1021_BANGS.fits.gz", "5866_BANGS.fits.gz")
    #for file in os.listdir(results_dir):
    #    if file.endswith("BANGS.fits.gz"):
    #        file_list.append(file)

    name = os.path.join(BeagleDirectories.results_dir, args.summary_config) 
    my_photometry.summary_catalogue.compute(file_list, name)

    # *****************************************************
    # *********** "BANGS MultiNest catalogue" ****************
    # *****************************************************

    file_name = "UVUDF_MultiNest.cat"

    # ********* Load ***************
    #my_photometry.multinest_catalogue.load(file_name)

    # ********* Compute ***************
    #file_list = list()
    #for file in os.listdir(results_dir):
    #    if file.endswith("MNstats.dat"):
    #        file_list.append(file)

    #n_par = 6
    #my_photometry.multinest_catalogue.compute( n_par, file_list, file_name)

    # *****************************************************
    # *********** Posterior Predictive Checks  ****************
    # *****************************************************

    file_name = "PPC.fits"

    # ********* Load ***************
    #my_photometry.PPC.load( file_name) 

    # ********* Compute ***************
    my_photometry.PPC.compute(my_photometry.observed_catalogue, 
            my_photometry.filters, 
            file_name=file_name)

    # ********* Plot ***************

    #my_photometry.PPC.plot_chi2()

    my_photometry.PPC.plot_p_value(broken_axis=True)

    stop

    # *****************************************************
    # *********** Residual Photometry  ****************
    # *****************************************************

    my_photometry.residual.compute(my_photometry.observed_catalogue,
            my_photometry.summary_catalogue, my_photometry.filters)

    # *****************************************************
    # *********** PDF  ****************
    # *****************************************************
