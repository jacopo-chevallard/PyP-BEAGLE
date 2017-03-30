#!/usr/bin/env python
import sys
import os
import argparse
import re
import numpy as np
import ConfigParser
import logging
from matplotlib import rc
from astropy.io import ascii
from astropy.io import fits

from pyp_beagle import *
from pyp_beagle import __version__
from pathos.multiprocessing import ProcessingPool 


def main():

    # Load the default command line argument parser
    parser = standard_parser()

    # Add package version
    parser.add_argument('-v', '--version', 
            action='version', 
            version='%(prog)s ' + __version__
            )

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

    # Check if the parameter file contains a LIST OF SPECTRA
    has_spectra = config.has_option('main', 'LIST OF SPECTRA')

    # Check if the parameter file contains a PHOTOMETRIC CATALOGUE
    has_photometry = config.has_option('main', 'PHOTOMETRIC CATALOGUE')

    # Global font size
    font = {'size': 16}
    rc('font', **font)

    # Get list of results files and object IDs from the results directory
    file_list, IDs = get_files_list()
    regex = re.compile(r"_MC\w+", re.IGNORECASE)

    # Load mock catalogue
    mock_catalogue = None
    if args.mock_file_name is not None:

        # JSON file containing the configuration for the mock catalogue plots
        params_file = os.path.join(BeagleDirectories.results_dir, args.json_file_mock)
        mock_catalogue = BeagleMockCatalogue(params_file, ignore_string=regex)
        mock_catalogue.load(args.mock_file_name)

    # JSON file containing the parameters to be plotted in the triangle plot
    params_file = os.path.join(BeagleDirectories.results_dir, args.json_file_triangle)

    # Set parameter names and labels
    my_PDF = PDF(params_file, 
            mock_catalogue=mock_catalogue)

    # Compute the summary catalogue
    if args.compute_summary:
        summary_catalogue = BeagleSummaryCatalogue()
        summary_catalogue.compute(file_list)

    # Comparison plots of true vs retrieved values 
    if args.mock_file_name is not None:
        if not summary_catalogue.exists():
            summary_catalogue.compute(file_list)
        summary_catalogue.load()

        mock_catalogue.compare_hist(summary_catalogue, overwrite=True)
        mock_catalogue.compare(summary_catalogue, overwrite=True)

    # ---------------------------------------------------------
    # --------- Post-processing of photometric data -----------
    # ---------------------------------------------------------
    if has_photometry:

        # Initialize an instance of the main "Photometry" class
        my_photometry = Photometry(key=args.ID_key, x_log=args.plot_log_wl)

        # We can load a set of photometric filters
        filters_file = os.path.expandvars(config.get('main', 'FILTERS FILE'))
        my_photometry.filters.load(filters_file)

        # Load observed photometric catalogue
        file_name = os.path.expandvars(config.get('main', 'PHOTOMETRIC CATALOGUE'))
        my_photometry.observed_catalogue.load(file_name)

    # ---------------------------------------------------------
    # -------- Post-processing of spectroscopic data ----------
    # ---------------------------------------------------------
    if has_spectra:

        # Initialize an instance of the main "Spectrum" class
        my_spectrum = Spectrum(params_file, 
                resolution=args.resolution, 
                plot_line_labels=args.plot_line_labels, 
                mock_catalogue=mock_catalogue)

        my_spectrum.observed_spectrum.configure(config=config)


        # File containing list of input spectra
        inputSpectraFileName = os.path.expandvars(config.get('main', 'LIST OF SPECTRA'))
        inputSpectraFile = open(inputSpectraFileName, 'r')

        file_names = list()

        for ID in IDs:

            # Plot the "triangle plot"
            #print "ID: ", ID

            for line in inputSpectraFile:
                # Get rid of the "\n" char at the end of the line
                line = line.strip()
                line = os.path.join(os.path.dirname(inputSpectraFileName), line)
                if ID in line:
                    file_names.append(line)
                    break

    # Create "pool" of processes
    if args.np > 1:
        pool = ProcessingPool(nodes=args.np)

    # Plot the marginal SED
    if args.plot_marginal:
        if args.np > 1:
            if has_spectra:
                pool.map(my_spectrum.plot_marginal, IDs, file_names)

            if has_photometry:
                pool.map(my_photometry.plot_marginal, IDs)
        else:
            for i, ID in enumerate(IDs):
                if has_spectra:
                    my_spectrum.plot_marginal(ID, file_names[i])

                if has_photometry:
                    my_photometry.plot_marginal(ID)

    # Plot the triangle plot
    if args.plot_triangle:
        if args.np > 1:
            pool.map(my_PDF.plot_triangle, IDs)
        else:
            for ID in IDs:
                my_PDF.plot_triangle(ID)
