# -*- coding: utf-8 -*-
"""
Created on Fri Jan 27 08:38:07 2023

@author: Niko Tutland
"""

import os
qf_options = {}
OG_PATH = os.getcwd()
PROJ_PATH = os.path.join(OG_PATH,"7.QUICFIRE-MODEL","projects")

###Settings to change###

# Simulation parameters
qf_options['QFVD'] = 5 #version 4 or 5
PROJ_FOLDER = "LandisTester" #folder containing QF inputs
qf_options['SimTime'] = 900
qf_options['print_times'] = 30

# Domain settings
qf_options['nx'] = 200 #starting values, may be modified in code
qf_options['ny'] = 200 #"
qf_options['nz'] = 22  #"

# Constant surface fuel parameters
qf_options['fuel_moisture'] = 0.10
qf_options['fuel_height'] = 0.5

# Topo settings
qf_options['topo_custom'] = False
qf_options['max_topo'] = 0 #set for flat right now, may use custom topo in the future

# Wind settings
qf_options['windspeed'] = 2.235 #m/s
qf_options['winddir'] = 270

# Ignitions settings
qf_options['custom_ig'] = True
qf_options['ig_xmin'] = 395.0 #m (Adam's original pattern here. Modified in code for LANDIS runs)
qf_options['ig_ymin'] = 2.0
qf_options['ig_xlen'] = 5.0
qf_options['ig_ylen'] = 397.0
# Driptorch options (for custom_ig = True)
qf_options['ig_method'] = "total" #options: "aerial", "drip", "total". Determines number of energy packets at ignitions
qf_options['ig_pattern'] = "strip" #options: "strip", "flank", "head", "back", "ring"
qf_options['ig_firingdir'] = 270
qf_options['ig_rate'] = 1.5 #m/s
qf_options['ig_fast'] = False #override rate and for fast (unrealistic) ignitions?
qf_options['ig_linetype'] = "line" #options: "line", "dash", "dot"
qf_options['ig_dashlen'] = None #length of ignition line (m) only for dash type
qf_options['ig_gaplen'] = None #length of gaps between ignitions (m) for dash and dot. When linetype = dashed, None means gap will be same as dash
qf_options['ig_depth'] = 60 #space between strip or flank igniters (m)
qf_options['ig_crewsize'] = 8 #crew size for strip or flank
qf_options['ig_offset'] = 20 #distance to firing area boundary for head or back (m)

###DO NOT CHANGE###
qf_options['RUN_PATH'] = os.path.join(PROJ_PATH,PROJ_FOLDER)
qf_options['dx'] = qf_options['dy'] = 2 #m
qf_options['dz'] = 1