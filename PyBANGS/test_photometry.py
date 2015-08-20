import os

from bangs_photometry import Photometry
from bangs_utils import BangsDirectories

results_dir = "/Users/efcl/BANGS/results/ECL13/scratch"
BangsDirectories.results_dir = results_dir

my_photometry = Photometry()

my_photometry.filters.load(os.path.expandvars("$BANGS_FILTERS/filters_ECL13_colName.dat"))


file_name = os.path.expandvars("$BANGS_DATA/ECL13/ECL13.cat")
my_photometry.observed_catalogue.load(file_name)

print my_photometry.observed_catalogue.data[my_photometry.observed_catalogue.data['ID'] == '18.0']

file_list=[]
for file in os.listdir(results_dir):
    if file.endswith("BANGS.fits.gz"):
        file_list.append(file)

#my_photometry.summary_catalogue.compute(file_list, file_name)

#my_photometry.residual.compute(my_photometry.observed_catalogue, 
#        my_photometry.summary_catalogue, my_photometry.filters)

my_photometry.plot_marginal('18.0')
