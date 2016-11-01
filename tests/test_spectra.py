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
import beagle_multiprocess

from pathos.multiprocessing import ProcessingPool 


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
    name = os.path.join(BeagleDirectories.results_dir, 'BEAGLE-input-files', param_file) 
    config.read(name)

    # File containing list of input spectra
    inputSpectraFileName = os.path.expandvars(config.get('main', 'LIST OF SPECTRA'))
    inputSpectraFile = open(inputSpectraFileName, 'r')

    # Global font size
    font = {'size': 16}
    rc('font', **font)

    # Get list of results files and object IDs from the results directory
    file_list, IDs = get_files_list()

    params_file = os.path.join(BeagleDirectories.results_dir, "params_names.json")
    # Initialize an instance of the main "Spectrum" class
    my_spectrum = Spectrum(params_file)

    my_spectrum.observed_spectrum.configure(config=config)

    # Set parameter names and labels
    my_PDF = PDF(params_file)

    #file_name = "BEAGLE_summary_catalogue.fits"
    #name = os.path.join(BeagleDirectories.results_dir, args.summary_config) 
    #my_spectrum.summary_catalogue.compute(file_list, name)

    my_spectrum.summary_catalogue.load()
    my_spectrum.mock_catalogue.load()
    my_spectrum.mock_catalogue.compare(my_spectrum.summary_catalogue, overwrite=True)

    stop

    observations_list = list()
    file_names = list()
    mock_file_names = list()
    ID_names = list()

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
                ID_names.append(ID)
                mock_file_name = line.split("_BEAGLE_")[0] + "_BEAGLE_MAP.fits.gz"
                mock_file_names.append(mock_file_name)
                break


    my_spectrum.mock_catalogue.compute(mock_file_names, overwrite=True)
    stop

    #p=Pool(processes=1)
    print "ID_names: ", ID_names[0]
    print "file_names: ", file_names[0]
    pool = ProcessingPool(nodes=16)
    kwargs = {'observation_name':file_names}
    args = ID_names
    #pool.map(my_spectrum.plot_marginal, args, kwargs)
    #pool.map(my_spectrum.plot_marginal, ID_names, observation_name=file_names)
    pool.map(my_PDF.plot_triangle, ID_names, mock_file_names)
    #stop

    #for ID, line in observations_list:
       # args, kwargs = ID, {'observation_name':line}
       # p.apply_async(my_spectrum.plot_marginal, args, kwargs)
        #my_spectrum.plot_marginal(ID, observation_name=line)
        #mock_file_name = line.split("_BEAGLE_")[0] + "_BEAGLE_MAP.fits.gz"
        #my_PDF.plot_triangle(ID, mock_file_name=mock_file_name)
