import os
import ast
import numpy as np


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

class PhotometricFilters:

    def __init__(self):

        self.index = list()
        self.colName = list()
        self.label = list()
        self.errcolName = list()

    def load(self, fileName):

        # read the filter file used to run BANGS

        for line in open(fileName):
            if line.strip() and not line.startswith("#"):
                # Split the line into different columns using whitespaces
                if 'units:value:' in line:
                    self.units = ast.literal_eval(line.split('units:value:')[1])
                elif 'units:' in line:
                    self.units = jy_to_erg(line.split('units:')[1])
                elif 'index:' in line:
                    self.index.append(line.split('index:')[1].split()[0])
                    self.colName.append(line.split('flux:colName:')[1].split()[0])
                    self.errcolName.append(line.split('fluxerr:colName:')[1].split()[0])
                    self.label.append(line.split('label:')[1].split()[0])

        self.n_bands = len(self.index)

        # Read filter transmission functions
        filt_log = os.path.expandvars("$BANGS_FILTERS/filters.log")
        n_wl_points = list()
        for line in open(filt_log):
            n_wl_points.append(int(line.split()[-1]))

        n_wl_points = np.array(n_wl_points)

        self.start_line = list()
        self.end_line = list()
        for i, indx in enumerate(self.index):
            self.start_line.append(int(np.sum(n_wl_points[0:int(indx)-1]) + int(indx)) )
            self.end_line.append(int(self.start_line[i]+n_wl_points[int(indx)-1]))

        filt_trans = os.path.expandvars("$BANGS_FILTERS/filterfrm.res")
        lines = open(filt_trans).readlines()

        self.wl_eff = list()

        for i in range(len(self.index)):
            n_wl = self.end_line[i]-self.start_line[i]
            wl = np.zeros(n_wl)
            t_wl = np.zeros(n_wl)
            for j, l in enumerate(lines[self.start_line[i]+1:self.end_line[i]]):
                wl[j], t_wl[j] = l.split()

            self.wl_eff.append(np.sum(wl*t_wl) / np.sum(t_wl))

