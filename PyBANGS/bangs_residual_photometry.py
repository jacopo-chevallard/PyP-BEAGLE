import cPickle

class ResidualPhotometry:

    def load_residuals(self, filename):
        """ 
        Load the pre-computed residul photometry from a cPickle file.

        Parameters
        ----------
        filename : str
            Name of the file containing the cPicke dump.
        """ 

        try:
            file = open(filename, 'rb')
            self.residual_kde_pdf = cPickle.load(file)
            file.close()
        except IOError as e:
            print "I/O error({0}): {1}".format(e.errno, e.strerror)

    def compute_residuals(self, observed_catalogue, bangs_summary_catalogue,
            filters, summary_stat=None, cPickleName=None):

        """ 
        Compute the residual photometry between an observed and model catalogues.  

        Parameters
        ----------
        observed_catalogue : Table object
            

        bangs_summary_catalogue : Table object
            

        filters : PhotometricFilters class object
            
        summary_stat : str, optional
            Either 'mean' or 'median', determines which type of summary
            statistics is used to compute the residual. By default the median
            of the marginal PDF is used.

        """

        _summary_stat = "median"
        if summary_stat is not None:
            _summary_stat = summary_stat

        bangs_data = bangs_summary_catalogue.hdulist['MARGINAL PHOTOMETRY'].data

        catalogue_data = observed_catalogue.data

        # As a sanity check, check if the ID match among the two catalogues, by
        # random picking some indices here and there...
        if (bangs_data['ID'][0] != bangs_data['ID'][0]) or
        (bangs_data['ID'][-1] != bangs_data['ID'][-1]):

        self.residual_kde_pdf = list() 

        jy = 1.E-26

        for i in range(filters.n_bands):

            obs_flux = catalogue_data[filters.colName[i]] * filters.units / jy
            obs_flux_err = catalogue_data[filters.errcolName[i]] * filters.units  / jy

            name = '_' + filters.label[i] + '_'
            model_flux = bangs_data[name+'_'+_summary_stat] / jy
            model_flux_err = 0.5 * (bangs_data[name+'_68.00'][:,1]-bangs_data[name+'_68.00'][:,0]) / jy

            mask = np.zeros(len(obs_flux), dtype=bool)
            mask[(obs_flux > 0.) & (obs_flux_err > 0.) & (model_flux > 0.) & (model_flux_err > 0.)] = True

            # Compute the residual in magnitude scale
            residual = -2.5*np.log10(model_flux[mask]/obs_flux[mask])

            # The residual error follows from standard error propagation
            residualErr = np.sqrt((2.5/model_flux[mask]*model_flux_err[mask])**2 + (2.5/obs_flux[mask]*obs_flux_err[mask])**2)

            # Compute the distribution of the residual by means of a weighted kernel density estimation
            kde_pdf = WeightedKDE.gaussian_kde(residual, weights = 1./residualErr)

            self.residual_kde_pdf.append(kde_pdf)

        if cPickleName is not None:
            cPickle.dump(self.residual_kde_pdf, cPickleName, cPickle.HIGHEST_PROTOCOL)

    def plot(self):

        # Compute interpolant of cumulative PDF (linear interpolation)
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        x_grid = np.linspace(-5.0,5.0,5000.)
        kde_pdf_grid = np.array(kde_pdf(x_grid))

        cumul_pdf = np.cumsum(kde_pdf_grid)
        # Be sure the it sums up to 1
        cumul_pdf /= cumul_pdf[-1]

        # Compute interpolant of cumulative PDF (linear interpolation)
        interp_cumul_pdf = interp1d(cumul_pdf, x_grid)

        print "filter: ", filters.label[i]
        print 'median:', interp_cumul_pdf(0.5)

        #print "x_grid: ", x_grid
        #print "y_val: ", y_val
        ax.plot(x_grid, kde_pdf_grid , ls="-", color = "black")
        plt.show()
        plt.close(fig)


