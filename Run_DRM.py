# -*- coding: utf-8 -*-
"""
Created on Wed Jul  5 09:03:41 2023

@author: Niko Tutland

Script for running the DRM for a L2-QF deomnstration paper. This assumes each
LANDIS run is already completed. 
"""

import sys,os
import os.path
import subprocess
import shutil
from time import sleep
import numpy as np
from QUICFire_options import qf_options
import DRM_framework_coupling as DRM
sys.path.insert(0, "1.LANDIS-MODEL")
import Driptorch_Support as dts
import LANDIS_to_Treelist as Landis
import Run_LANDIS as Run
import Crop_LANDIS as Crop


def main(nyears,ncycyear,ncycle,OG_PATH,VDM):
    build_trees(OG_PATH)
    cycle = 0      # current iteration (will be looped through range(0,ncycle))
    # Build Landis Parameters object for spinup
    L2_params = Run.LandisParams(OG_PATH, nyears, ncycyear, ncycle, cycle, spinup=True)
    os.chdir("..")
    # Crop to fire domain
    epsg = Crop.Landis(L2_params)
    Treelist_params = Landis.toTreelist(L2_params)
    # Define fire grid size and ignitions
    build_ignitions(qf_options,Treelist_params,OG_PATH,epsg)
    nfuel = Treelist_params.num_spp
    # Run Trees and QUIC-Fire
    qf_options['nx'], qf_options['ny'], qf_options['nz'] = Treelist_params.nx, Treelist_params.ny, Treelist_params.nz
    DRM.print_qf_inputs(qf_options)
    DRM.runTreeQF(VDM,FM,nfuel,qf_options['nx'],qf_options['ny'],qf_options['nz'],1)      # runs the tree program to create QF inputs
    ft_bragg_ff(OG_PATH,divide_by=16)
    DRM.runQF(0,VDM,qf_options)
    # Get fire effects on trees
    LiveDead = np.array(DRM.runCrownScorch(1,VDM,L2_params,Treelist_params))
    LiveDead = np.insert(LiveDead,0,1)

def build_ignitions(qf_options,Treelist_params,OG_PATH,epsg):
    if qf_options['custom_ig']:
        dts.write_ignitions(shape_path = os.path.join(OG_PATH,"1.LANDIS-MODEL","Shapefiles","burn_plot.shp"), 
                            ig_path = os.path.join(qf_options['RUN_PATH']), 
                            firing_dir = qf_options['ig_firingdir'], 
                            burn_path = os.path.join(OG_PATH,"1.LANDIS-MODEL","Shapefiles","burn_domain.shp"), 
                            epsg = epsg,
                            rate = qf_options['ig_rate'], 
                            fast = qf_options['ig_fast'], 
                            line_type = qf_options['ig_linetype'], 
                            dash_length = qf_options['ig_dashlen'], 
                            gap_length = qf_options['ig_gaplen'],
                            pattern = qf_options['ig_pattern'], 
                            crew_size = qf_options['ig_crewsize'], 
                            depth = qf_options['ig_depth'], 
                            offset = qf_options['ig_offset'])
    else:
        qf_options['ig_xmin'] = 100
        qf_options['ig_ymin'] = Treelist_params.ny/2 # half of half the y length
        qf_options['ig_xlen'] = 10
        qf_options['ig_ylen'] = Treelist_params.ny #since dy is 2, this is half the length of the y side of the domain
    print("ignitions written to QF run")

def build_trees(OG_PATH):
    os.chdir("5.TREES-QUICFIRE")
    if os.path.exists(os.path.join(OG_PATH,"5.TREES-QUICFIRE","treelist_VDM.dat")):
        with subprocess.Popen(
            ["wsl","make","clean"], stdout=subprocess.PIPE
        ) as process:

            def poll_and_read():
                print(f"{process.stdout.read1().decode('utf-8')}")
            
            while process.poll() != 0:
                poll_and_read()
                sleep(1)
            if process.poll()==0:
                print('make clean successful - running make')
                with subprocess.Popen(
                    ["wsl","make"], stdout=subprocess.PIPE
                ) as process:
                
                    def poll_and_read():
                        print(f"{process.stdout.read1().decode('utf-8')}")
                    
                    while process.poll() != 0:
                        poll_and_read()
                        sleep(1)
                    if process.poll()==0:
                        print('trees successfully compiled')
    else:
        with subprocess.Popen(
            ["wsl","make"], stdout=subprocess.PIPE
        ) as process:
        
            def poll_and_read():
                print(f"{process.stdout.read1().decode('utf-8')}")
            
            while process.poll() != 0:
                poll_and_read()
                sleep(1)
            if process.poll()==0:
                print('trees successfully compiled')
    os.chdir("..")
    
def ft_bragg_ff(OG_PATH,divide_by):
    fpath = os.path.join(OG_PATH,"1.LANDIS-MODEL/VDM2FM","VDM_litter_trees.dat")
    shutil.copyfile(src=fpath,dst=os.path.join(OG_PATH,"1.LANDIS-MODEL/VDM2FM","VDM_litter_trees_OG.dat"))
    rhof = np.loadtxt(fpath)
    rhof = np.divide(rhof,divide_by)
    np.savetxt(fpath,rhof)
    print("**WARNING: Surface fuel values divided by ", divide_by, "**")

if __name__=="__main__":
    nyears=50      # number of years for spinup and transient runs
    ncycyear=0    # number of cyclical year run (zero indicates no subsequent VDM runs)
    ncycle=1      # number of loops
    VDM = "LANDIS"
    FM = "QUICFIRE"
    OG_PATH = os.getcwd()
    main(nyears,ncycyear,ncycle,OG_PATH,VDM)