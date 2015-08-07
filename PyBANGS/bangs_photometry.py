


from bangs_filters import PhotometricFilters


class ObservedCatalogue:

    def load_catalogue(self, filename):

        """ 
        Load a photometric catalogue of observed sources. It automatically
        detects, and loads, FITS or ASCII files depending on the suffix.

        Parameters
        ----------
        filename : str
            Contains the file name of the catalogue.

        """

        try:
            if filename.endswith(('fits', 'fit', 'FITS', 'FIT')):
                self.data = fits.open(filename)[1].data
                self.columns = fits.open(filename)[1].columns
            else:
                self.data = ascii.read(filename, Reader=ascii.basic.CommentedHeader)
        except IOError as e:
            print "I/O error({0}): {1}".format(e.errno, e.strerror)


class Photometry:

    def load_filters(self, filename):

        """ 
        Load various information about a set of photometric filters used
        during a BANGS run. 

        Parameters
        ----------
        filename : str
            Contains the filter file used in the BANGS run.

        Notes
        -----
        For consistency with BANGS (and to mnimize errors), it uses the
        '$BANGS_FILTERS' environment variable to load the 'filters.log' and
        'filterfrm.res' files.
        """

        self.filters = PhotometricFilters()
        self.filters.load_filters(filters_filename)

    def load_observed_catalogue(self, filename):

        """ 
        Load a photometric catalogue of observed sources. It automatically
        detects, and loads, FITS or ASCII files depending on the suffix.

        Parameters
        ----------
        filename : str
            Contains the file name of the catalogue.

        """

        self.observed_catalogue = ObservedCatalogue()
        self.observed_catalogue.load_catalogue(filename)

    def load_bangs_summary_catalogue(self, filename):

        self.bangs_summary_catalogue = BangsSummaryCatalogue()
        self.bangs_summary_catalogue.load_catalogue(filename)


