import sys
sys.path.append("../PyP-BEAGLE")
import os
import argparse
import ConfigParser
import logging
from matplotlib import rc
from astropy.io import ascii
from astropy.io import fits

from beagle_spectra import Spectrum
from beagle_pdf import PDF
from beagle_utils import BeagleDirectories, get_files_list, find_file
from beagle_parsers import standard_parser


if __name__ == '__main__':

    # 
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

    # File containing list of input spectra
    inputSpectraFileName = os.path.expandvars(config.get('main', 'LIST OF SPECTRA'))
    inputSpectraFile = open(inputSpectraFileName, 'r')

    # Global font size
    font = {'size': 16}
    rc('font', **font)

    # Get list of results files and object IDs from the results directory
    file_list, IDs = get_files_list()

    # Initialize an instance of the main "Spectrum" class
    my_spectrum = Spectrum()

    my_spectrum.observed_spectrum.configure(config=config)

    # Set parameter names and labels
    my_PDF = PDF(os.path.join(BeagleDirectories.results_dir, "params_names.json"))

    for ID in IDs:

        # Plot the "triangle plot"
        print "ID: ", ID

        for line in inputSpectraFile:
            # Get rid of the "\n" char at the end of the line
            line = line.strip()
            line = os.path.join(os.path.dirname(inputSpectraFileName), line)
            if ID in line:
                my_spectrum.observed_spectrum.load(line)
                my_spectrum.plot_marginal(ID)
                break

        mock_file_name = line.split("_BEAGLE_")[0] + "_BEAGLE_MAP.fits.gz"
        my_PDF.plot_triangle(ID, mock_file_name=mock_file_name)
