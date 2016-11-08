
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
        '-v', '--verbose',
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
        '-s', '--summary-config',
        help="JSON file containing the configuration for the computation of the summary catalogue",
        action="store", 
        type=str, 
        dest="summary_config",
        default="summary_config.json"
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
