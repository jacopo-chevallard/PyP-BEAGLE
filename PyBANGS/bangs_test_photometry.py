import os
import argparse
import logging
from astropy.io import ascii
from astropy.io import fits
from bangs_photometry import Photometry

# 
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
    action="store_const", dest="loglevel", const=logging.INFO,
)
args = parser.parse_args()    
logging.basicConfig(level=args.loglevel)

# Initialize an instance of the main "Photometry" class
results_dir = "/Users/jchevall/My_Codes/BANGS_REPOS/BANGS_root/results/UVUDF/RAND_0.2/mass_sfr_cb14"
my_photometry = Photometry(results_dir)

# We now have access to other classes!

# *****************************************************
# *********** Filters ****************
# *****************************************************

# We can load a set of photometric filters
my_photometry.filters.load("/Users/jchevall/My_Codes/BANGS_REPOS/BANGS_root/filters/filters_UVUDF.dat")

# *****************************************************
# *********** Observed Catalogue ****************
# *****************************************************

# ********** Loading *****************
file_name = "/Users/jchevall/My_Codes/BANGS_REPOS/BANGS_root/data/UVUDF/hlsp_uvudf_hst_v2.0_cat_Types_no_star.fits"
my_photometry.observed_catalogue.load(file_name)

# *****************************************************
# *********** "BANGS summary catalogue" ****************
# *****************************************************

file_name = os.path.join(results_dir, "BANGS_summary_catalogue_TEST.fits")

# ********* Test loading of the catalogue ***************
#my_photometry.summary_catalogue.load(file_name)

# ********* Test computation of the catalogue ***************
file_list = list()
for file in os.listdir(results_dir):
    if file.endswith("BANGS.fits.gz"):
        file_list.append(file)

my_photometry.summary_catalogue.compute(file_list, results_dir, file_name)

# *****************************************************
# *********** "BANGS MultiNest catalogue" ****************
# *****************************************************

# ********* Test computation of the catalogue ***************
file_list = list()
for file in os.listdir(results_dir):
    if file.endswith("MNstats.dat"):
        file_list.append(file)

file_name = "UVUDF_MultiNest.cat"
n_par = 6
my_photometry.multinest_catalogue.compute( n_par, file_list, results_dir, file_name)

# *****************************************************
# *********** Posterior Predictive Checks  ****************
# *****************************************************

# then you can compute the Posterior Predictive Checks
file_name = "PPC.fits"
my_photometry.PPC.compute(my_photometry.observed_catalogue, 
        my_photometry.filters, 
        my_photometry.results_dir, 
        file_name=file_name)

# *****************************************************
# *********** Residual Photometry  ****************
# *****************************************************

# *****************************************************
# *********** PDF  ****************
# *****************************************************
