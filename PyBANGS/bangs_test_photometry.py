import os
import argparse
import logging
from astropy.io import ascii
from astropy.io import fits

from bangs_photometry import Photometry
from bangs_pdf import PDF
from bangs_utils import BangsDirectories

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
results_dir = "/Users/jchevall/Coding/BANGS/files/root/results/UVUDF/RAND_0.2/mass_specific_sfr_cb14"
BangsDirectories.results_dir = results_dir

my_photometry = Photometry()
#ID = 5866

#my_PDF = PDF(os.path.join(results_dir, "params_names.json"))
#my_PDF.plot_triangle(ID)

# We now have access to other classes!

# *****************************************************
# *********** Filters ****************
# *****************************************************

# We can load a set of photometric filters
my_photometry.filters.load(os.path.expandvars("$BANGS_FILTERS/filters_UVUDF.dat"))

# *****************************************************
# *********** Observed Catalogue ****************
# *****************************************************

# ********** Loading *****************
#file_name = os.path.expandvars("$BANGS_DATA/UVUDF/hlsp_uvudf_hst_v2.0_cat_Types_no_star.fits")
#my_photometry.observed_catalogue.load(file_name)


# ********** Plotting of the marginal photometry *****************
#my_photometry.plot_marginal(ID)

# *****************************************************
# *********** "BANGS summary catalogue" ****************
# *****************************************************

file_name = "BANGS_summary_catalogue_TEST.fits"

# ********* Load ***************
#my_photometry.summary_catalogue.load(file_name)

# ********* Compute ***************
file_list = ("1021_BANGS.fits.gz", "5866_BANGS.fits.gz")
#for file in os.listdir(results_dir):
#    if file.endswith("BANGS.fits.gz"):
#        file_list.append(file)

my_photometry.summary_catalogue.compute(file_list, file_name)

stop

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

#my_photometry.PPC.plot_p_value()

# *****************************************************
# *********** Residual Photometry  ****************
# *****************************************************

my_photometry.residual.compute(my_photometry.observed_catalogue,
        my_photometry.summary_catalogue, my_photometry.filters)

# *****************************************************
# *********** PDF  ****************
# *****************************************************
