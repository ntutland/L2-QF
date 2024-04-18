# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 09:11:58 2023

@author: Niko Tutland
"""

import os
import sys
import numpy as np
import drawfire
from read_inputs import read_fireca_field

def main(prj_folder):
    calc_rhof_change(prj_folder)

def calc_rhof_change(prj_folder):
    qu, qf, ignitions, flags, fb, output_folder = drawfire.import_inputs(prj_folder)
    fuel_dens = read_fireca_field("fuels-dens-", qf.ntimes, qf.time, qf, 0, prj_folder)
    with np.errstate(divide = 'ignore'):
        Drhofp = np.divide(fuel_dens[-1],fuel_dens[0])
    Drhofp[np.isnan(Drhofp)] = 1
    PercRhofChange = 0.0
    print("Assembling PercentFuelChange.txt")
    fname10 = "PercentFuelChange.txt"
    with open(fname10, 'w') as pfc:
        for locz in range(qf.nz):
            for locy in range(qf.ny):
                for locx in range(qf.nx):
                    PercRhofChange = Drhofp[locy,locx,locz]
                    pfc.write('{}\n'.format(PercRhofChange))
    del fuel_dens

if __name__ == '__main__':
    prj_folder = os.getcwd() if len(sys.argv == 1) else sys.argv[1]   
    main(prj_folder)