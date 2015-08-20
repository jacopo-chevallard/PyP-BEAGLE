import os
from bangs_utils import BangsDirectories
import bangs_pdf

results_dir = "/Users/efcl/BANGS/results/ECL13/scratch"
BangsDirectories.results_dir = results_dir

pdf = bangs_pdf.PDF("/Users/efcl/BANGS/results/ECL13/scratch/params_names.json")

#pdf.plot_triangle('18.0')

os.chdir(results_dir)

fileList = os.listdir(results_dir)

for files in fileList:
    if '.fits.gz' in files:
        print files
        id = files.replace('_BANGS.fits.gz', '')
        print files, id
        if id != '2.0' and id != '10.0':
            pdf.plot_triangle(id)
