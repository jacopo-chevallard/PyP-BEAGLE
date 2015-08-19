import os
from bangs_utils import BangsDirectories
import bangs_pdf

results_dir = "/Users/efcl/BANGS/results/ECL13/scratch"
BangsDirectories.results_dir = results_dir

pdf = bangs_pdf.PDF("/Users/efcl/BANGS/results/ECL13/first_TEST/params.json")

pdf.plot_triangle('18.0')
