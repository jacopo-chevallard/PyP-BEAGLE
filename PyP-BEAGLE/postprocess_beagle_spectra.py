#!/usr/bin/env python
import sys
import os
sys.path.append(os.path.join(os.environ['PYP_BEAGLE'], "PyP-BEAGLE"))
import argparse
import re
import numpy as np
import ConfigParser
import logging
from matplotlib import rc
from astropy.io import ascii
from astropy.io import fits

from beagle_mock_catalogue import BeagleMockCatalogue
from beagle_spectra import Spectrum
from beagle_pdf import PDF
from beagle_utils import BeagleDirectories, get_files_list, find_file
from beagle_parsers import standard_parser
import beagle_multiprocess

from pathos.multiprocessing import ProcessingPool 


if __name__ == '__main__':

    # Load the default command line argument parser
    parser = standard_parser()

    # Get parsed arguments
    args = parser.parse_args()    

    # Set logging level
    logging.basicConfig(level=args.loglevel)

    # Set directory containing BEAGLE results files
    BeagleDirectories.results_dir = args.results_dir

    # Read parameter file
    config = ConfigParser.SafeConfigParser()

    # Check if you passed a name of the parameter file, otherwise search for a
    # suitable parameter file in the BEAGLE-input-files folder
    param_file = None
    if args.param_file is not None:
        param_file = args.param_file
    else:
        for file in os.listdir(os.path.join(BeagleDirectories.results_dir, 
            BeagleDirectories.beagle_input_files)):
            if file.endswith('.param'):
                file_name = os.path.join(BeagleDirectories.results_dir, 
                        BeagleDirectories.beagle_input_files, file)
                if '[main]' in open(file_name).readline():
                    param_file = file
                    break

    if param_file is None:
        raise ValueError("No parameter file has been specified, nor a parameter file can be " 
                "found in the Beagle results directory")

    # Set name of parameter file
    BeagleDirectories.param_file = param_file 

    # Parse the parameter file
    name = os.path.join(BeagleDirectories.results_dir, BeagleDirectories.beagle_input_files, 
            param_file)
    config.read(name)

    # File containing list of input spectra
    inputSpectraFileName = os.path.expandvars(config.get('main', 'LIST OF SPECTRA'))
    inputSpectraFile = open(inputSpectraFileName, 'r')

    # Global font size
    font = {'size': 16}
    rc('font', **font)

    # Get list of results files and object IDs from the results directory
    file_list, IDs = get_files_list()
    regex = re.compile(r"_MC\w+", re.IGNORECASE)
    #test = ('V_150_MC_0_snr_PS_CL', 'V_10_MC_0_snr_PS_CLE')
    #for t in test:
    #    tt = regex.sub('', t)
    #    print "tt: ", tt
    # Load mock catalogue
    mock_catalogue = None
    if args.mock_file_name is not None:

        # JSON file containing the configuration for the mock catalogue plots
        params_file = os.path.join(BeagleDirectories.results_dir, args.json_file_mock)
        mock_catalogue = BeagleMockCatalogue(params_file, ignore_string=regex)
        mock_catalogue.load(args.mock_file_name)

    # JSON file containing the parameters to be plotted in the triangle plot
    params_file = os.path.join(BeagleDirectories.results_dir, args.json_file_triangle)

    # Initialize an instance of the main "Spectrum" class
    my_spectrum = Spectrum(params_file, 
            resolution=100, 
            plot_line_labels=True, 
            mock_catalogue=mock_catalogue)

    my_spectrum.observed_spectrum.configure(config=config)

    # Set parameter names and labels
    my_PDF = PDF(params_file, 
            mock_catalogue=mock_catalogue)

    # Compute the summary catalogue
    if args.compute_summary:
        my_spectrum.summary_catalogue.compute(file_list)

    # Comparison plots of true vs retrieved values 
    if args.mock_file_name is not None:
        if not my_spectrum.summary_catalogue.exists():
            my_spectrum.summary_catalogue.compute(file_list)
        my_spectrum.summary_catalogue.load()

        
        class_indices = list() ; class_color = list()
        z_min = 4 ; z_max = 8 ; 
        redshift = mock_catalogue.hdulist['GALAXY PROPERTIES'].data['redshift']
        indx = np.where((redshift < 6.))[0]
        class_indices.append(indx)
        indx = np.where((redshift >= 6.))[0]
        class_indices.append(indx)

        mock_catalogue.compare_hist(my_spectrum.summary_catalogue, overwrite=True, class_indices=class_indices)
        mock_catalogue.compare(my_spectrum.summary_catalogue, overwrite=True)

    observations_list = list()
    file_names = list()

    for ID in IDs:

        # Plot the "triangle plot"
        #print "ID: ", ID

        for line in inputSpectraFile:
            # Get rid of the "\n" char at the end of the line
            line = line.strip()
            line = os.path.join(os.path.dirname(inputSpectraFileName), line)
            if ID in line:
                observations_list.append((ID, line))
                file_names.append(line)
                break

    # Create "pool" of processes
    if args.np > 1:
        pool = ProcessingPool(nodes=args.np)

    # Plot the marginal SED
    if args.plot_marginal:
        if args.np > 1:
            pool.map(my_spectrum.plot_marginal, IDs, file_names)
        else:
            for ID, file in zip(IDs, file_names):
                my_spectrum.plot_marginal(ID, file)

    # Plot the triangle plot
    if args.plot_triangle:
        if args.np > 1:
            pool.map(my_PDF.plot_triangle, IDs)
        else:
            for ID in IDs:
                my_PDF.plot_triangle(ID)
