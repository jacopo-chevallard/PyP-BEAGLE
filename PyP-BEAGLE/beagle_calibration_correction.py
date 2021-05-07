import json
from collections import OrderedDict
import numpy as np
import scipy.special as special
    
class CalibrationCorrection:

    def __init__(self):
    
        self.type=None
        self.degree=None
        self.has_correction=False
        
    def configure(self, config, params_file):


        self.has_correction = config.has_option('main', 'FLUX CALIBRATION CORRECTION')
        
        if self.has_correction and params_file is not None:
        
            #Names of parameters and whether they are fitted or not, supplying the
            #values if they are fixed
            with open(params_file) as f:
                self.coeff_params = json.load(f, object_pairs_hook=OrderedDict)
            line = config.get('main', 'FLUX CALIBRATION CORRECTION')
            if 'type:' in line:
                self.type = line.split("type:")[1].split()[0]
                #Type can currently only be Polynomial
                if self.type != "polynomial":
                    raise ValueError("Calibration Correction of type `" + type + "` not recognised!")
            if 'degree:' in line:
                self.degree = np.int(line.split("degree:")[1].split()[0])
            
    
    def return_correction(self, x, coeff):
    
        y = np.zeros_like(x)
        if self.type == "polynomial":
            for i in range(self.degree+1):
                y = y + coeff[i]*np.power(x,i)
        if self.type == "legendre":
            for i in range(self.degree+1):
                y = y + coeff[i]*special.eval_legendre(i,x)
        return y

