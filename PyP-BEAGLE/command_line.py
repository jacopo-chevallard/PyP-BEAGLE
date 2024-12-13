#!/usr/bin/env python
from __future__ import absolute_import
import sys
import os
import argparse
import re
import numpy as np
import six.moves.configparser
import logging
from matplotlib import rc
from astropy.io import ascii
from astropy.io import fits

from .beagle_parsers import standard_parser
from .beagle_utils import (
    BeagleDirectories,
    get_files_list,
    configure_matplotlib,
    trimFitsSuffix,
)
from .beagle_mock_catalogue import BeagleMockCatalogue
from .beagle_summary_catalogue import BeagleSummaryCatalogue
from .beagle_photometry import Photometry
from .beagle_filters import PhotometricFilters
from .beagle_spectral_indices import SpectralIndices
from .beagle_spectra import Spectrum
from .beagle_pdf import PDF

from ._version import __version__
from six.moves import zip

# Only need to set spawn method on macOS (Darwin)
if sys.platform == "darwin":
    from pathos.helpers import mp

    mp.set_start_method("spawn", force=True)

from pathos.multiprocessing import ProcessingPool


# Configure matplotlib for the main process
configure_matplotlib()


# Create a wrapper function that configures matplotlib for each worker
def worker_init():
    configure_matplotlib()


def main():

    # Load the default command line argument parser
    parser = standard_parser()

    # Add package version
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="%(prog)s " + __version__ + " - Author: Jacopo Chevallard",
    )

    # Get parsed arguments
    args = parser.parse_args()

    # Convert argparse.ArgumentParser() Namespace object into a dictionary (see https://stackoverflow.com/a/16878364)
    args_dict = vars(args)

    # Set logging level
    logging.basicConfig(level=args.loglevel)

    # Set directory containing BEAGLE results files
    BeagleDirectories.results_dir = args.results_dir

    # Configure matplotlib
    configure_matplotlib()

    # Read parameter file
    config = six.moves.configparser.SafeConfigParser(strict=False)

    # Set fontsize
    BeagleDirectories.fontsize = args.fontsize
    BeagleDirectories.inset_fontsize_fraction = args.inset_fontsize_fraction
    font = {"size": BeagleDirectories.fontsize}
    rc("font", **font)

    # Check if you passed a name of the parameter file, otherwise search for a
    # suitable parameter file in the BEAGLE-input-files folder
    param_file = None
    if args.param_file is not None:
        param_file = args.param_file
    else:
        for file in os.listdir(
            os.path.join(
                BeagleDirectories.results_dir, BeagleDirectories.beagle_input_files
            )
        ):
            if file.endswith(".param"):
                file_name = os.path.join(
                    BeagleDirectories.results_dir,
                    BeagleDirectories.beagle_input_files,
                    file,
                )
                if "[main]" in open(file_name).readline():
                    param_file = file
                    break

    if param_file is None:
        raise ValueError(
            "No parameter file has been specified, nor a parameter file can be "
            "found in the Beagle results directory"
        )

    # Set name of parameter file
    BeagleDirectories.param_file = param_file

    # Parse the parameter file
    name = os.path.join(
        BeagleDirectories.results_dir, BeagleDirectories.beagle_input_files, param_file
    )
    config.read(name)

    # Check if the parameter file contains a LIST OF SPECTRA
    has_spectra = config.has_option("main", "LIST OF SPECTRA")

    # Check if the parameter file contains a PHOTOMETRIC CATALOGUE
    has_photometry = config.has_option("main", "PHOTOMETRIC CATALOGUE")

    # Check if the parameter file contains a SPECTRAL INDICES CATALOGUE
    has_spec_indices = config.has_option("main", "SPECTRAL INDICES CATALOGUE")

    # Get list of results files and object IDs from the results directory
    file_list, IDs = get_files_list(suffix=args.suffix)
    if len(file_list) == 0:
        raise ValueError(
            "No Beagle results files are present in the directory "
            + BeagleDirectories.results_dir
        )

    regex = None
    if args.regex_ignore is not None:
        # regex = re.compile(r"_MC\w+", re.IGNORECASE)
        regex = re.compile(args.regex_ignore, re.IGNORECASE)

    for i, (ID, file) in enumerate(zip(IDs, file_list)):
        if args.ID_list is not None:
            if not ID in args.ID_list:
                IDs.remove(ID)
                file_list.remove(file)
                continue
        if regex is not None:
            _ID = regex.sub("", ID)
            IDs[i] = _ID

    # Load mock catalogue
    mock_catalogue = None
    if args.mock_file_name is not None:

        # JSON file containing the configuration for the mock catalogue plots
        params_file = os.path.join(BeagleDirectories.results_dir, args.json_file_mock)
        mock_catalogue = BeagleMockCatalogue(
            params_file, ignore_string=regex, plot_title=args.plot_title
        )
        mock_catalogue.load(args.mock_file_name)

    # JSON file containing the parameters to be plotted in the triangle plot
    params_file = os.path.join(BeagleDirectories.results_dir, args.json_file_triangle)

    # Compute the summary catalogue
    if args.compute_summary:

        summary_catalogue = BeagleSummaryCatalogue(
            credible_intervals=args.credible_interval,
            config_file=args.summary_config,
            n_proc=args.n_proc,
            overwrite=args.overwrite,
            flatten_columns=args.flatten_columns,
        )

        if not summary_catalogue.exists() or args.overwrite:
            summary_catalogue.compute(file_list)
        summary_catalogue.load()

    if args.latex_table_params is not None:
        if not args.compute_summary:
            parser.error(
                "The `--latex-table-params` option requires the `--compute-summary` option to be set"
            )
        summary_catalogue.make_latex_table(
            args.latex_table_params,
            IDs=args.ID_list,
            average_errors=args.average_errors,
            column_wise=args.latex_column,
        )

    if args.extract_MAP:
        if not args.compute_summary:
            parser.error(
                "The `--extract-MAP` option requires the `--compute-summary` option to be set"
            )
        summary_catalogue.extract_MAP_solution(file_list)

    # Comparison plots of true vs retrieved values
    if args.mock_file_name is not None:
        # mock_catalogue.plot_input_param_distribution()
        mock_catalogue.compare_hist(summary_catalogue, overwrite=True)
        mock_catalogue.compare(summary_catalogue, overwrite=True)

    # ---------------------------------------------------------
    # --------- Post-processing of photometric data -----------
    # ---------------------------------------------------------
    if has_photometry and args.plot_marginal:

        # We can load a set of photometric filters
        try:
            filters_file = os.path.expandvars(config.get("main", "FILTERS FILE"))
            filters_throughputs = None
        except:
            filters_file = os.path.expandvars(
                config.get("main", "FILTERS CONFIGURATION")
            )
            try:
                filters_throughputs = os.path.expandvars(
                    config.get("main", "FILTERS THROUGHPUTS")
                )
            except:
                filters_throughputs = None

        my_filters = PhotometricFilters()

        my_filters.load(
            filters_file,
            filters_folder=args.filters_folder,
            filters_throughputs=filters_throughputs,
        )

        # Initialize an instance of the main "Photometry" class
        my_photometry = Photometry(my_filters, **args_dict)

        # Load observed photometric catalogue
        file_name = os.path.expandvars(config.get("main", "PHOTOMETRIC CATALOGUE"))
        my_photometry.observed_catalogue.load(file_name)

    # ---------------------------------------------------------
    # --------- Post-processing of spectral indices data -----------
    # ---------------------------------------------------------
    if has_spec_indices:

        # Initialize an instance of the main "SpectralIndices" class
        my_spec_indices = SpectralIndices(**args_dict)

        file_name = os.path.expandvars(
            config.get("main", "SPECTRAL INDICES CONFIGURATION")
        )
        my_spec_indices.observed_catalogue.configure(file_name)

        # Load observed catalogue
        file_name = os.path.expandvars(config.get("main", "SPECTRAL INDICES CATALOGUE"))
        my_spec_indices.observed_catalogue.load(file_name)

    # ---------------------------------------------------------
    # -------- Post-processing of spectroscopic data ----------
    # ---------------------------------------------------------
    if has_spectra and args.plot_marginal:

        # Initialize an instance of the main "Spectrum" class
        my_spectrum = Spectrum(**args_dict)

        my_spectrum.observed_spectrum.configure(config=config)

        # File containing list of input spectra
        inputSpectraFileName = os.path.expandvars(config.get("main", "LIST OF SPECTRA"))

        lines = list()
        with open(inputSpectraFileName, "r") as f:
            for line in f:
                # Get rid of the "\n" char at the end of the line
                line = line.strip()
                line = os.path.join(os.path.dirname(inputSpectraFileName), line)
                lines.append(line)

        spec_file_names = dict()

        for ID in IDs:
            _ID = ID
            if regex is not None:
                _ID = regex.sub("", _ID)
            for line in lines:
                _line = trimFitsSuffix(os.path.basename(line))
                if regex is not None:
                    _line = regex.sub("", _line)
                _line = _line.split("_BEAGLE")[0]
                if _ID == _line:
                    spec_file_names[ID] = line
                    break

    # Create "pool" of processes
    if args.n_proc > 1:
        pool = ProcessingPool(nodes=args.n_proc, initializer=worker_init)

    # Plot the marginal SED
    if args.plot_marginal:
        if args.n_proc > 1:
            if has_spectra:
                pool.map(
                    my_spectrum.plot_marginal,
                    [ID for ID in IDs if ID in spec_file_names],
                    [spec_file_names[ID] for ID in IDs if ID in spec_file_names],
                )

            if has_photometry:
                pool.map(my_photometry.plot_marginal, IDs)

            if has_spec_indices:
                pool.map(my_spec_indices.plot_line_fluxes, IDs)
        else:
            for i, ID in enumerate(IDs):
                if has_spectra and ID in spec_file_names:
                    my_spectrum.plot_marginal(ID, spec_file_names[ID])

                if has_photometry:
                    my_photometry.plot_marginal(ID)

                if has_spec_indices:
                    my_spec_indices.plot_line_fluxes(ID)

    # Plot the triangle plot
    if args.plot_triangle:

        # Set parameter names and labels
        my_PDF = PDF(params_file, mock_catalogue=mock_catalogue, **args_dict)

        if args.n_proc > 1:
            pool.map(my_PDF.plot_triangle, IDs)
        else:
            for ID in IDs:
                my_PDF.plot_triangle(ID)
