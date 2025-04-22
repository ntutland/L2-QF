# -*- coding: utf-8 -*-
"""
Created on Fri Jan 27 08:38:07 2023

@author: Niko Tutland
"""

import os

qf_options = {}
OG_PATH = os.getcwd()
PROJ_PATH = os.path.join(OG_PATH, "7.QUICFIRE-MODEL", "projects")

###Settings to change###

# Simulation parameters
qf_options["RUN_SIM"] = False
qf_options["QFVD"] = "v6"  # "v5" or "v6"
PROJ_FOLDER = "LandisTester"  # folder containing QF inputs
qf_options["SimTime"] = 10800
qf_options["print_times"] = 30

# Domain settings
qf_options["nx"] = 200  # starting values, may be modified in code
qf_options["ny"] = 200  # "
qf_options["nz"] = 22  # "

# Constant surface fuel parameters
qf_options["fuel_moisture"] = 0.10
qf_options["fuel_height"] = 0.2

# Topo settings
qf_options["topo_custom"] = False
qf_options["max_topo"] = 0  # set for flat right now, may use custom topo in the future

# Wind settings
qf_options["windspeed"] = 2.23  # m/s
qf_options["winddir"] = 270

# Ignitions settings
qf_options["custom_ig"] = False
qf_options["ig_xmin"] = (
    395.0  # m (Adam's original pattern here. Modified in code for LANDIS runs)
)
qf_options["ig_ymin"] = 2.0
qf_options["ig_xlen"] = 5.0
qf_options["ig_ylen"] = 397.0
qf_options["ig_method"] = (
    "drip"  # options: "aerial", "drip", "total". Determines number of energy packets at ignitions
)
qf_options["fireline_width"] = 8

# Driptorch ignition options
qf_options["ig_firingdir"] = qf_options["winddir"]
qf_options["ig_rate"] = 1.5
qf_options["ig_fast"] = False
qf_options["ig_linetype"] = "dash"
qf_options["ig_dashlen"] = 4
qf_options["ig_gaplen"] = 4
qf_options["ig_pattern"] = "strip"  # options: "strip", "flank", "head", "back"
qf_options["ig_crewsize"] = 4
qf_options["ig_depth"] = 20
qf_options["ig_offset"] = 10


###DO NOT CHANGE###
qf_options["PROJ_PATH"] = PROJ_PATH
qf_options["dx"] = qf_options["dy"] = 2  # m
qf_options["dz"] = 1
