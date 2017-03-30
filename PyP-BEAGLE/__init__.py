from beagle_utils import *
from beagle_parsers import standard_parser
from beagle_filters import PhotometricFilters
from beagle_photometry import Photometry
from beagle_pdf import PDF
from beagle_spectra import Spectrum
from beagle_summary_catalogue import BeagleSummaryCatalogue
from beagle_mock_catalogue import BeagleMockCatalogue
from beagle_residual_photometry import ResidualPhotometry
from beagle_multinest_catalogue import MultiNestCatalogue
from beagle_posterior_predictive_checks import PosteriorPredictiveChecks

from pkg_resources import get_distribution, DistributionNotFound
import os.path

try:
    _dist = get_distribution('pyp_beagle')
    # Normalize case for Windows systems
    dist_loc = os.path.normcase(_dist.location)
    here = os.path.normcase(__file__)
    if not here.startswith(os.path.join(dist_loc, 'pyp_beagle')):
        # not installed, but there is another version that *is*
        raise DistributionNotFound
except DistributionNotFound:
    __version__ = 'Please install this project with setup.py'
else:
    __version__ = _dist.version
