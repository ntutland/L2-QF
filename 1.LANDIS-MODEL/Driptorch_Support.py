# -*- coding: utf-8 -*-
"""
Created on Fri Dec  9 12:25:43 2022

@author: ntutland
"""

import driptorch as dt
import os
import geopandas as gpd
import pickle


def write_ignitions(
    shape_path,
    ig_path,
    firing_dir,
    burn_path,
    epsg=4326,
    rate=1.5,
    fast=False,
    line_type="line",
    dash_length=None,
    gap_length=None,
    pattern="strip",
    crew_size=8,
    depth=60,
    offset=20,
):
    burn_plot = gpd.read_file(shape_path)
    burn_domain = gpd.read_file(burn_path)
    firing_area = setup_burnunit(burn_plot, firing_dir, epsg)
    burn_domain = setup_burndomain(burn_domain, firing_dir, epsg=epsg)
    igniter = create_igniters(rate, fast, line_type, dash_length, gap_length)
    crew = setup_crew(igniter, crew_size)
    ig_pattern = setup_ignitions(firing_area, crew, pattern, depth, offset)
    ig_pattern.to_quicfire(
        burn_domain,
        filename=os.path.join(ig_path, "ignite.dat"),
        dst_epsg=epsg,
        resolution=2,
    )
    delete_ignite_header(ig_path)
    if line_type == "dot":
        delete_ignite_zeros(ig_path)
    ig_pickle = open(os.path.join(ig_path, "ignite.obj"), "wb")
    pickle.dump(ig_pattern, ig_pickle)
    ig_pickle.close()
    how_long(ig_pattern)


def setup_burnunit(burn_plot, firing_dir, epsg):
    burn_unit = dt.BurnUnit(burn_plot.iloc[0, -1], firing_dir, utm_epsg=epsg)
    firing_area = burn_unit.buffer_control_line(10)
    return firing_area


def setup_burndomain(burn_domain, firing_dir, epsg):
    burn_domain = dt.BurnUnit(burn_domain.iloc[0, -1], firing_dir, utm_epsg=epsg)
    return burn_domain


def create_igniters(rate, fast, line_type, dash_length, gap_length):
    fastrate = 10  # m/s
    if fast:
        rate = fastrate
    if line_type == "line":
        igniter = dt.Igniter(rate)
    elif line_type == "dash":
        igniter = dt.Igniter(rate, dash_length, gap_length)
    elif line_type == "dot":
        igniter = dt.Igniter(rate, gap_length=gap_length)
    else:
        print("Unsupported ignition line type")

    return igniter


def setup_crew(igniter, size):
    crew = dt.IgnitionCrew.clone_igniter(igniter, size)
    return crew


def setup_ignitions(burn_unit, crew, pattern, depth=60, offset=20):
    if pattern == "strip":
        ignitions = dt.firing.Strip(burn_unit, crew)
        ig_pattern = ignitions.generate_pattern(spacing=0, depth=depth)
    elif pattern == "flank":
        ignitions = dt.firing.Flank(burn_unit, crew)
        ig_pattern = ignitions.generate_pattern(depth=depth)
    elif pattern == "head":
        ignitions = dt.firing.Head(burn_unit, crew)
        ig_pattern = ignitions.generate_pattern(offset=offset)
    elif pattern == "back":
        ignitions = dt.firing.Back(burn_unit, crew)
        ig_pattern = ignitions.generate_pattern(offset=offset)
    elif pattern == "ring":
        ignitions = dt.firing.Ring(burn_unit, crew)
        ig_pattern = ignitions.generate_pattern(offset=offset)
    else:
        print("Unsupported ignition pattern")

    return ig_pattern


def delete_ignite_header(path):
    file = os.path.join(path, "ignite.dat")
    # list to store file lines
    lines = []
    # read file
    with open(file, "r") as fp:
        # read an store all lines into list
        lines = fp.readlines()
    # Write file
    with open(file, "w") as fp:
        # iterate each line
        for number, line in enumerate(lines):
            if number not in range(7):
                fp.write(line)


def delete_ignite_zeros(path):
    file = os.path.join(path, "ignite.dat")
    # list to store file lines
    lines = []
    # read file
    with open(file, "r") as fp:
        # read an store all lines into list
        lines = fp.readlines()
        linecount = len(lines) - 6  # first 6 lines aren't ignitions

    deleted = 0
    # Write file
    with open(file, "w") as fp:
        # iterate each line
        for line in lines:
            if line.startswith("0"):
                deleted += 1
            else:
                fp.write(line)

    with open(file, "r") as fp:
        # read an store all lines into list
        lines = fp.readlines()

    # Write file
    with open(file, "w") as fp:
        # iterate each line
        for line in lines:
            if line.find(" 0 ") != -1:
                deleted += 1
            else:
                fp.write(line)

    newlinecount = linecount - deleted
    with open(file, "r") as fp:
        # read an store all lines into list
        lines = fp.readlines()
    with open(file, "w") as fp:
        for line in lines:
            if line.startswith("naerial="):
                fp.write("naerial={}\n".format(newlinecount))
            else:
                fp.write(line)


def how_long(ignitions):
    t = ignitions.elapsed_time
    print("Ignitions take {} seconds".format(t))


if __name__ == "__main__":
    OG_PATH = os.getcwd()
    shape_path = os.path.join(OG_PATH, "Shapefiles", "burn_plot.shp")
    ig_path = os.path.join(OG_PATH, "Ignitions")
    firing_dir = 270

    write_ignitions(
        shape_path, ig_path, firing_dir, line_type="dash", dash_freq=10, pattern="head"
    )
