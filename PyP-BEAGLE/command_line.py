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

from pathos.multiprocessing import ProcessingPool 
from pyp_beagle import *
from _version import __version__


def main():

    # Load the default command line argument parser
    parser = standard_parser()

    # Add package version
    parser.add_argument('-v', '--version', 
            action='version', 
            version='%(prog)s ' + __version__ + ' - Author: Jacopo Chevallard'
            )

    # Get parsed arguments
    args = parser.parse_args()    

    # Set logging level
    logging.basicConfig(level=args.loglevel)

    # Set directory containing BEAGLE results files
    BeagleDirectories.results_dir = args.results_dir

    # Read parameter file
    config = ConfigParser.SafeConfigParser()

    # Set fontsize
    BeagleDirectories.fontsize = args.fontsize
    BeagleDirectories.inset_fontsize_fraction = args.inset_fontsize_fraction
    font = {'size': BeagleDirectories.fontsize}
    rc('font', **font)

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

    # Get list of results files and object IDs from the results directory
    file_list, IDs = get_files_list(suffix=args.suffix)
    if len(file_list) == 0:
        raise ValueError("No Beagle results files are present in the directory " + BeagleDirectories.results_dir)

    regex = None
    if args.regex_ignore is not None:
        #regex = re.compile(r"_MC\w+", re.IGNORECASE)
        regex = re.compile(args.regex_ignore, re.IGNORECASE)

    if args.ID_list is not None:
        IDs_ = IDs
        file_list_ = file_list
        for ID, file in zip(IDs_, file_list_):
            if regex is not None:
                ID_ = regex.sub('', ID)
            else:
                ID_ = ID
            if not ID_ in args.ID_list:
                IDs.remove(ID)
                file_list.remove(file)

    # Load mock catalogue
    mock_catalogue = None
    if args.mock_file_name is not None:

        # JSON file containing the configuration for the mock catalogue plots
        params_file = os.path.join(BeagleDirectories.results_dir, args.json_file_mock)
        mock_catalogue = BeagleMockCatalogue(params_file, ignore_string=regex, plot_title=args.plot_title)
        mock_catalogue.load(args.mock_file_name)

    # JSON file containing the parameters to be plotted in the triangle plot
    params_file = os.path.join(BeagleDirectories.results_dir, args.json_file_triangle)

    # Compute the summary catalogue
    summary_catalogue = BeagleSummaryCatalogue(credible_intervals=args.credible_interval, n_proc=args.n_proc)
    if args.compute_summary:
        if not summary_catalogue.exists():
            summary_catalogue.compute(file_list)

    if args.latex_table_params is not None:
        if not summary_catalogue.exists():
            summary_catalogue.compute(file_list)
        summary_catalogue.load()
        summary_catalogue.make_latex_table(args.latex_table_params, IDs=args.ID_list)

    if args.extract_MAP:
        summary_catalogue.extract_MAP_solution(file_list)

    # Comparison plots of true vs retrieved values 
    if args.mock_file_name is not None:
        if not summary_catalogue.exists():
            summary_catalogue.compute(file_list)
        summary_catalogue.load()

        #mock_catalogue.plot_input_param_distribution()
        mock_catalogue.compare_hist(summary_catalogue, overwrite=True)
        mock_catalogue.compare(summary_catalogue, overwrite=True)

    # ---------------------------------------------------------
    # --------- Post-processing of photometric data -----------
    # ---------------------------------------------------------
    if has_photometry:

        # Initialize an instance of the main "Photometry" class
        my_photometry = Photometry(key=args.ID_key, 
                x_log=args.plot_log_wl, 
                log_flux=args.plot_log_flux,
                plot_single_solution=args.plot_single_solution,
                plot_full_SED=args.plot_full_SED)

        # We can load a set of photometric filters
        try:
            filters_file = os.path.expandvars(config.get('main', 'FILTERS FILE'))
            filters_throughputs = None
        except:
            filters_file = os.path.expandvars(config.get('main', 'FILTERS CONFIGURATION'))
            filters_throughputs = os.path.expandvars(config.get('main', 'FILTERS THROUGHPUTS'))

        my_photometry.filters.load(filters_file, 
                filters_folder=args.filters_folder, 
                filters_throughputs=filters_throughputs)

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
                plot_full_SED=args.plot_full_SED,
                wl_range=args.wl_range,
                line_labels_json=args.line_labels_json,
                plot_line_labels=args.plot_line_labels, 
                mock_catalogue=mock_catalogue,
                log_flux=args.plot_log_flux,
                show_residual=args.show_residual,
                print_ID=args.print_ID,
                wl_rest=args.wl_rest,
                draw_steps=args.draw_steps,
                plot_suffix=args.plot_suffix,
                wl_units=args.wl_units)

        my_spectrum.observed_spectrum.configure(config=config)


        # File containing list of input spectra
        inputSpectraFileName = os.path.expandvars(config.get('main', 'LIST OF SPECTRA'))

        lines = list()
        with open(inputSpectraFileName, 'r') as f:
            for line in f:
                # Get rid of the "\n" char at the end of the line
                line = line.strip()
                line = os.path.join(os.path.dirname(inputSpectraFileName), line)
                lines.append(line)

        file_names = list()

        for ID in IDs:
            ID_ = ID
            if regex is not None:
                ID_ = regex.sub('', ID_)
            for line in lines:
                line_ = trimFitsSuffix(os.path.basename(line))
                if regex is not None:
                    line_ = regex.sub('', line_)
                if ID_ == line_:
                    file_names.append(line)
                    break

    # Create "pool" of processes
    if args.n_proc > 1:
        pool = ProcessingPool(nodes=args.n_proc)

    # Plot the marginal SED
    if args.plot_marginal:
        if args.n_proc > 1:
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

        # Set parameter names and labels
        my_PDF = PDF(params_file, 
                mock_catalogue=mock_catalogue,
                plot_single_solution=args.plot_single_solution)

        if args.n_proc > 1:
            pool.map(my_PDF.plot_triangle, IDs)
        else:
            for ID in IDs:
                my_PDF.plot_triangle(ID)
