# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 08:04:05 2022

@author: Niko Tutland
"""

opdict = {}

# Run info
opdict["states"] = [
    "NC",
    "VA",
    "TN",
    "SC",
    "GA",
    "FL",
]  # states containing and bordering study area
opdict["fia_spec"] = [
    "PITA",
    "PIPA2",
    "PIEC2",
    "PIEL",
    "QUAL",
    "QULA2",
    "LIST2",
    "ACRU",
    "LITU",
]  # fia species symbols for tree species in study area
opdict["landis_spec"] = [
    "LobPine",
    "LongleafPine",
    "ShortPine",
    "SlashPine",
    "WhiteOak",  # landis designation for tree species in study area
    "TurkeyOak",
    "SweetGum",
    "RedMaple",
    "TulipTree",
]  # (input in same order as fia_spec)
opdict["region_flag"] = (
    3  # where is the AOI? 1 = California, 2 = other western state, 3 = midwest, eastern, or southern state
)
opdict["age_bin"] = 10
opdict["aoi_elev"] = (
    325  # elevation of the study area in feet [should we calculate this in the code?]
)

# Fuels
opdict["bulk_density"] = 0.7  # constant canopy bulk density (kg/m^3)
opdict["cl_factor"] = 0.8  # percent of crown above point of maximum crown diameter (m)
opdict["moisture"] = 1  # constant canopy moisture content
opdict["sizescale"] = 0.0005  # canopy fuel radius (nominal fuel size - solids) (m)

# Spatial info
opdict["crop_domain"] = True

# Fire effects
opdict["mortality_thresholds"] = [
    0.7,
    0.3,
    0.7,
    0.7,
    0.8,  # threshold for percent fuel remaining after fire sim, under which the tree does not survive
    0.8,
    0.9,
    0.9,
    0.9,
]  # (input in same order as fia_spec and landis_spec)

## Make sure the community output biomass extension prints outputs at an interval that will capture both the spinup end year and the cycles' end year
