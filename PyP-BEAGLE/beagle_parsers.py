
import argparse
import logging

def standard_parser():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        '-d', '--debug',
        help="Print lots of debugging statements",
        action="store_const", 
        dest="loglevel", 
        const=logging.DEBUG,
        default=logging.WARNING
    )

    parser.add_argument(
        '--verbose',
        help="Be verbose",
        action="store_const", 
        dest="loglevel", 
        const=logging.INFO
    )

    parser.add_argument(
        '-r', '--results-dir',
        help="Directory containing BEAGLE results",
        action="store", 
        type=str, 
        dest="results_dir", 
        required=True
    )

    parser.add_argument(
        '-p', '--parameter-file',
        help="Parametr file used in the BEAGLE run",
        action="store", 
        type=str, 
        dest="param_file"
    )

    parser.add_argument(
        '--json-triangle',
        help="JSON file used for the triangle plots.",
        action="store", 
        type=str, 
        dest="json_file_triangle",
        default="params_names.json"
    )

    parser.add_argument(
        '--plot-marginal',
        help="Plot the marginal SED.",
        action="store_true", 
        dest="plot_marginal" 
        )

    parser.add_argument(
        '--plot-single-solution',
        help="FITS table containing the columns 'ID' and 'row_index' identifying the solution, "\
                "among those available in the posterior PDF, to overplot on the triangle and marginal plots.",
        action="store", 
        type=str, 
        dest="plot_single_solution" 
        )

    parser.add_argument(
        '--log-wavelength',
        help="Plot logarithmic axis for wavelength.",
        action="store_true", 
        dest="plot_log_wl" 
        )

    parser.add_argument(
        '--plot-line-labels',
        help="Label the strongest emission lines.",
        action="store_true", 
        dest="plot_line_labels" 
        )

    parser.add_argument(
        '--plot-title',
        help="Add title to plot.",
        action="store", 
        type=str, 
        dest="plot_title" 
        )

    parser.add_argument(
        '--plot-full-SED',
        help="Overplot the full resolution model SED.",
        action="store_true", 
        dest="plot_full_SED" 
        )

    parser.add_argument(
        '--wl-range',
        help="Wavelength range(s) (in micron) to plot.",
        action="store",
        type=float,
        nargs='+',
        dest="wl_range" 
        )

    parser.add_argument(
        '--spectral-resolution',
        help="(Approximate) resolution of the spectra, used to determine which emission lines \
                to label (lines which are blended are not labelled separately).",
        action="store",
        type=float,
        dest="resolution" 
        )

    parser.add_argument(
        '--plot-triangle',
        help="Plot the triangle plot.",
        action="store_true", 
        dest="plot_triangle" 
        )

    parser.add_argument(
        '--compute-summary',
        help="Compute the summary catalogue.",
        action="store_true", 
        dest="compute_summary" 
        )

    parser.add_argument(
        '--extract-MAP',
        help="Extract the Maximum-a-Posteriori solution.",
        action="store_true", 
        dest="extract_MAP" 
        )

    parser.add_argument(
        '--json-summary',
        help="JSON file containing the configuration for the computation of the summary catalogue",
        action="store", 
        type=str, 
        dest="summary_config",
        default="summary_config.json"
    )

    parser.add_argument(
        '--mock-catalogue',
        help="Mock catalogue containing the 'true' parameter values.",
        action="store", 
        type=str, 
        dest="mock_file_name"
    )

    parser.add_argument(
        '--json-mock',
        help="JSON file used for the mock true vs retrieved parameters plots.",
        action="store", 
        type=str, 
        dest="json_file_mock",
        default="params_names.json"
    )

    parser.add_argument(
        '--ID-key',
        help="Name of the column containing the object ID in the source catalogue",
        action="store", 
        type=str, 
        dest="ID_key",
        default="ID"
    )

    # Number of processors to use in the multi-processor parts of the analysis
    parser.add_argument(
        '-np',
        help="Number of processors to use",
        action="store", 
        type=int, 
        dest="np",
        default=1
    )

    return parser
