import os
import logging
import ast
import numpy as np
from astropy.table import Table, Column
from astropy.io import fits, ascii


class UnitsError(Exception):

    def __init__(self, value):
        self.value = value
        self.msg = "Invalid unit: ", self.value

def jy_to_erg(jy):

    if jy.strip().lower() == 'jy':
        return 1.e-23
    elif jy.strip().lower() == 'millijy':
        return 1.e-23*1.e-03
    elif jy.strip().lower() == 'microjy':
        return 1.e-23*1.e-06
    elif jy.strip().lower() == 'nanojy':
        return 1.e-23*1.e-09
    else:
        try:
            # Raise an exception with argument
            raise UnitsError(jy)
        except UnitsError, arg:
            # Catch the custom exception
            print 'Error: ', arg.msg

class PhotometricFilters(object):

    def load(self, file_name, 
            filters_folder="$BEAGLE_FILTERS", 
            filters_throughputs=None):
        """ 
        Load various information about a set of photometric filters used
        during a BEAGLE run. 

        Parameters
        ----------
        file_name : str
            Contains the filter file used in the BEAGLE run.

        Notes
        -----
        For consistency with BEAGLE (and to mnimize errors), it uses the
        '$BEAGLE_FILTERS' environment variable to load the 'filters.log' and
        'filterfrm.res' files.
        """

        logging.info("Loading the `PhotometricFilters` file: " + file_name)

        # Count number of bands defined in filter file
        self.n_bands = 0
        with open(file_name) as f:
            for line in f:
                if ('index:' in line or 'name:' in line) and not line.startswith("#"):
                    self.n_bands += 1

        transmission = list()
        index = Column(name='index', dtype=np.int32, length=self.n_bands)
        name = Column(name='name', dtype='S40', length=self.n_bands)
        fileName = Column(name='fileName', dtype='S250', length=self.n_bands)
        wl_eff = Column(name='wl_eff', dtype=np.float32, length=self.n_bands)
        colName = Column(name='flux_colName', dtype='S20', length=self.n_bands)
        errcolName = Column(name='flux_errcolName', dtype='S20', length=self.n_bands)
        label = Column(name='label', dtype='S40', length=self.n_bands)

        min_rel_err = Column(name='min_rel_err', dtype=np.float32, length=self.n_bands)

        old_API = False

        # read the filter file used to run BEAGLE
        i = 0
        for line in open(file_name):
            if line.strip() and not line.startswith("#"):
                # Split the line into different columns using whitespaces
                if 'units:value:' in line:
                    self.units = ast.literal_eval(line.split('units:value:')[1])
                elif 'flux:conversion:' in line:
                    self.units = ast.literal_eval(line.split('flux:conversion:')[1])
                elif 'units:' in line:
                    self.units = jy_to_erg(line.split('units:')[1])
                elif 'name:' in line or 'index:' in line:
                    if 'index:' in line:
                        old_API = True
                        index[i] = line.split('index:')[1].split()[0]
                    else:
                        name[i] = line.split('name:')[1].split()[0]
                    if 'fileName:' in line:
                        fileName[i] = line.split('fileName:')[1].split()[0]
                    if 'flux:colName:' in line:
                        colName[i] = line.split('flux:colName:')[1].split()[0]
                    if 'fluxerr:colName:' in line:
                        errcolName[i] = line.split('fluxerr:colName:')[1].split()[0]
                    if 'fluxErr:colName:' in line:
                        errcolName[i] = line.split('fluxErr:colName:')[1].split()[0]
                    if 'label:' in line:
                        label[i] = line.split('label:')[1].split()[0]
                    min_rel_err[i] = 0.
                    if 'min_rel_err:' in line:
                        min_rel_err[i] = line.split('min_rel_err:')[1].split()[0]
                    i += 1    

        if old_API:
            # Read filter transmission functions
            filt_log = os.path.expandvars(os.path.join(filters_folder, "filters.log"))
            n_wl_points = list()
            for line in open(filt_log):
                n_wl_points.append(int(line.split()[-1]))

            n_wl_points = np.array(n_wl_points)

            self.start_line = list()
            self.end_line = list()
            for i, indx in enumerate(index):
                self.start_line.append(int(np.sum(n_wl_points[0:int(indx)-1]) + int(indx)) )
                self.end_line.append(int(self.start_line[i]+n_wl_points[int(indx)-1]))

            filt_trans = os.path.expandvars(os.path.join(filters_folder, "filterfrm.res"))
            lines = open(filt_trans).readlines()

            for i in range(len(index)):
                n_wl = self.end_line[i]-self.start_line[i]
                wl = np.zeros(n_wl)
                t_wl = np.zeros(n_wl)
                for j, l in enumerate(lines[self.start_line[i]:self.end_line[i]-1]):
                    wl[j], t_wl[j] = l.split()

                wl_eff[i] = np.sum(wl*t_wl) / np.sum(t_wl)

                transmission.append({'wl':wl, 't_wl':t_wl})
        else:
            hdulist = fits.open(filters_throughputs)

            for i in range(self.n_bands):

                if fileName[i].strip():
                    data = ascii.read(fileName[i])
                    wl, t_wl = data.field(0), data.field(1)
                else:
                    d = hdulist['TRANSMISSION'].data[name[i]][0]
                    wl, t_wl = d[0,:], d[1,:]

                wl_eff[i] = np.sum(wl*t_wl) / np.sum(t_wl)
                transmission.append({'wl':wl, 't_wl':t_wl})
            hdulist.close()

        trans = Column(name='transmission', data=transmission)
        my_cols = [index, name, colName, errcolName, label, wl_eff, min_rel_err, trans]
        
        self.columns = my_cols
        self.data = Table(my_cols)
