
from astropy.io import ascii
from astropy.io import fits
from bangs_photometry import Photometry

# Initialize an instance of the main "Photometry" class
my_photometry = Photometry("/Users/jchevall/My_Codes/BANGS_REPOS/BANGS_root/results/UVUDF/RAND_0.2/mass_sfr_cb14")

# We now have access to other classes!

# We can load a set of photometric filters
my_photometry.filters.load("/Users/jchevall/My_Codes/BANGS_REPOS/BANGS_root/filters/filters_UVUDF.dat")

# and an observed photometric catalogue
file_name = "/Users/jchevall/My_Codes/BANGS_REPOS/BANGS_root/data/UVUDF/hlsp_uvudf_hst_v2.0_cat_Types_no_star.fits"
my_photometry.observed_catalogue.load(file_name)

# and a "BANGS summary catalogue"
file_name = "/Users/jchevall/My_Codes/BANGS_REPOS/BANGS_root/results/UVUDF/RAND_0.2/mass_specific_sfr_cb14/BANGS_summary_catalogue.fits"
my_photometry.summary_catalogue.load(file_name)

# and a "MultiNest" catalogue
file_name = "/Users/jchevall/My_Codes/BANGS_REPOS/BANGS_root/results/UVUDF/RAND_0.2/mass_specific_sfr_cb14/UVUDF_MultiNest.cat"
#def compute(self, n_par, file_list, file_name):
#my_photometry.multinest_catalogue.compute( file_name=file_name)

# then you can compute the Posterior Predictive Checks
results_dir = "/Users/jchevall/My_Codes/BANGS_REPOS/BANGS_root/results/UVUDF/RAND_0.2/mass_sfr_cb14"
my_photometry.PPC.compute(my_photometry.observed_catalogue, 
        my_photometry.filters, 
        my_photometry.results_dir, 
        file_name='TEST.fits')
