import os
import argparse
import logging
import ConfigParser
import multiprocessing as mp
import copy_reg
import types
from matplotlib import rc
from astropy.io import ascii
from astropy.io import fits
import numpy as np

from bangs_photometry import Photometry
from bangs_pdf import PDF
from bangs_utils import BangsDirectories, get_results_files, getattr_proxy
from bangs_analysis import wrap_make_photometry_plots, my_print_ID, my_other_print_ID

def _reduce_method(m):
    if m.im_self is None:
        return getattr, (m.im_class, m.im_func.func_name)
    else:
        return getattr, (m.im_self, m.im_func.func_name)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument(
        '-d', '--debug',
        help="Print lots of debugging statements",
        action="store_const", dest="loglevel", const=logging.DEBUG,
        default=logging.WARNING,
    )

    parser.add_argument(
        '-v', '--verbose',
        help="Be verbose",
        action="store_const", 
        dest="loglevel", 
        const=logging.INFO
    )

    # Number of processors to use in the multi-processor parts of the analysis
    parser.add_argument(
        '-r', '--results-dir',
        help="Directory containing BEAGLE results",
        action="store", 
        type=str, 
        dest="results_dir", 
        required=True
    )

    # Number of processors to use in the multi-processor parts of the analysis
    parser.add_argument(
        '-np',
        help="Number of processors to use",
        action="store", 
        type=int, 
        dest="np",
        default=1
    )

    # Get parsed arguments
    args = parser.parse_args()    

    # Set logging level
    logging.basicConfig(level=args.loglevel)

    # Set directory containing BEAGLE results files
    BangsDirectories.results_dir = args.results_dir

    # Read parameter file
    config = ConfigParser.SafeConfigParser()

    parameter_file = os.path.join(BangsDirectories.results_dir, 'UVUDF_BANGS_head.param')
    config.read(parameter_file)

    # Global font size
    font = {'size': 16}
    rc('font', **font)

    # Get list of results files and object IDs from the results directory
    file_list, IDs = get_results_files()

    # Initialize an instance of the main "Photometry" class
    my_photometry = Photometry()

    # Set parameter names and labels
    my_PDF = PDF(os.path.join(BangsDirectories.results_dir, "params_names.json"))

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
    my_photometry.observed_catalogue.load(file_name)

    #rand_IDs = np.random.choice(IDs, size=10)
    rand_IDs = (1021, 5866)

    for ID in rand_IDs:
        my_photometry.plot_marginal(ID, replot=True)
        my_PDF.plot_triangle(ID, M_star=True, replot=True)

    stop
    #IDs = (2930, 21930, 3413, 7113, 7930, 9806, 9972, 4930, 55, 930)

    #for ID in IDs:
    #    my_photometry.plot_marginal(ID)
    #    my_PDF.plot_triangle(ID)

    pool = mp.Pool(4)
    results = list()

    copy_reg.pickle(types.MethodType, _reduce_method)

    for ID in IDs[0:7]:
       # print "ID: ", ID
        #pool.apply_async(print_ID, args=[config])
        #pool.map(getattr_proxy(my_photometry, "print_ID"), ID)
        #pool.apply_async(getattr_proxy(my_photometry, "plot_marginal"), ID)
        #pool.apply_async(my_photometry.plot_marginal, args=[my_photometry.plot_marginal, ID])
        #pool.apply_async(wrap_make_photometry_plots, (ID, my_photometry, my_PDF))

        #pool.apply_async(my_print_ID, (ID,), dict(photometry=my_photometry, PDF=my_PDF))
        my_other_print_ID(ID, my_photometry)
        #pool.apply_async(my_other_print_ID, args=[ID, my_photometry], callback=results.append)

    while len(results) < 8:
        continue

    pool.close()


    pause


    # ********** Plotting of the marginal photometry *****************
    #my_photometry.plot_marginal(ID)
    my_photometry.plot_replicated_data(ID)
    stop

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

    #my_photometry.summary_catalogue.compute(file_list, file_name)

    # *****************************************************
    # *********** "BANGS MultiNest catalogue" ****************
    # *****************************************************

    file_name = "UVUDF_MultiNest.cat"

    # ********* Load ***************
    my_photometry.multinest_catalogue.load(file_name)

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
    my_photometry.PPC.load( file_name) 

    # ********* Compute ***************
    #my_photometry.PPC.compute(my_photometry.observed_catalogue, 
    #        my_photometry.filters, 
    #        file_name=file_name)

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
