from __future__ import absolute_import
from astropy.io import ascii
from astropy.io import fits

from .beagle_utils import is_FITS_file
        
class ObservedCatalogue(object):

    def load(self, file_name):

        """ 
        Load a catalogue of observed sources. It automatically
        detects, and loads, FITS or ASCII files depending on the suffix.

        Parameters
        ----------
        file_name : str
            Contains the file name of the catalogue.
        """

        if is_FITS_file(file_name):
            self.data = fits.open(file_name)[1].data
            self.columns = fits.open(file_name)[1].columns
        else:
            self.data = ascii.read(file_name, Reader=ascii.basic.CommentedHeader)


