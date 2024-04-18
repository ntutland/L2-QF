# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 12:37:26 2023

@author: FireScience
"""

import os
import sys
import copy
import math
import pylab
import matplotlib as mpl
import numpy as np
import class_def as cl
import drawfire
import gen_images as img
import misc
import read_inputs as read
from read_inputs import read_fireca_field

class PrintSpecs:
    def __init__(self, output_folder, img_specs, fuel_idx, no_fuel_idx, moist_planes, dens_planes):
        self.output_folder = output_folder
        self.img_specs = img_specs
        self.fuel_idx = fuel_idx
        self.no_fuel_idx = no_fuel_idx
        self.moist_planes = moist_planes
        self.dens_planes = dens_planes

def main(prj_folder):
    ## WHICH PLANES TO PLOT? ##
    dens_planes = [1] # planes are not zero-indexed
    moist_planes = [1]
    
    # read input files
    qu, qf, ignitions, flags, fb, output_folder = drawfire.import_inputs(prj_folder)
    
    ##----More flexible version of drawfire.plot_outputs----##
    
    gen_gif = 1 #always generage gif, but only for plane 1
    
    # Setting image specs
    img_specs = cl.ImgClass(gen_gif)
    img.set_image_specifications(img_specs)

    # Create folder to save images
    misc.create_plots_folder(gen_gif, img_specs, prj_folder, 'Plots')
    
    if flags.topo > 0:
        print("\t-terrain elevation")
        img.plot_terrain(qu, img_specs)
    
    if flags.isfire == 1:
        
        print("\t-fuel density field")
        fuel_idx, no_fuel_idx = drawfire.compress_fuel_dens(qf, flags, output_folder)
        
        print_specs = PrintSpecs(output_folder,img_specs, fuel_idx, no_fuel_idx, moist_planes, dens_planes)

        # ------- Firetech ignitions
        print("\t-initial ignitions")
        img.plot_ignitions(qf, fuel_idx, ignitions.hor_plane, qf.horizontal_extent, img_specs)

        # ------ Targets
        target(prj_folder, qf, print_specs)

        # ------- Fuel height (ground level only)
        fuel_height(qf, print_specs)

        # ------- Mass burnt (vertically-integrated)
        mass_burnt(flags, qf, print_specs)

        # ------- Fuel mass
        fuel_density(flags, qf, print_specs)

        # ------- Emissions
        emissions(flags, qf, print_specs)

        # ------- Radiation
        radiation(flags, qf, print_specs)

        # ------- Energy to atmosphere
        eng_to_atmos(flags, qf, print_specs)

        # -------  Fireca winds
        qf_winds(flags, qf, print_specs)

        # -------  QU winds (instantaneous)
        qu_winds_inst(flags, qf, qu, print_specs)

        # -------  QU winds (average)
        qu_winds_avg(flags, qf, qu, print_specs)

        # ------- w plumes
        w_plumes(flags, qf, qu, print_specs)

        # ------- Reaction rate
        reaction_rate(flags, qf, print_specs)
        
        # ------- Fuel moisture
        fuel_moisture(flags, qf, print_specs)

def target(prj_folder, qf, ps):
    targets = read.get_targets(prj_folder, ps.output_folder, qf)
    if targets.num > 0:
        print("\t-target distances")
        img.plot_targets(qf, targets, ps.img_specs, ps.output_folder, ps.fuel_idx, ps.no_fuel_idx)

def fuel_height(qf,ps):
    print("\t-ground level fuel height")
    ground_level_fuel_height = read.read_ground_fuel_height(qf, ps.output_folder)
    img.plot_fuelheight(qf, ground_level_fuel_height, ps.img_specs)
    del ground_level_fuel_height

def mass_burnt(flags,qf,ps):
    if flags.perc_mass_burnt == 1:
        print("\t-% mass burnt")
        perc_mass_burnt = read_fireca_field("mburnt_integ-", qf.ntimes, qf.time, qf, 1, ps.output_folder)
        img.plot_percmassburnt(qf, perc_mass_burnt, ps.no_fuel_idx, ps.img_specs, flags)
        del perc_mass_burnt
    return

def fuel_density(flags, qf, ps):
        if flags.fuel_density == 1:
            fuel_dens = read_fireca_field("fuels-dens-", qf.ntimes, qf.time, qf, 0, ps.output_folder)
            if qf.dx == 2:
                all_planes = ps.dens_planes
            else:
                all_planes = list([1])
            for iplane in all_planes:
                if iplane <= qf.nz:
                    print("\t-fuel mass, plane: %d" % iplane)
                    plot_2d_field(False, qf, iplane, fuel_dens, "Fuel density [kg/m^3]", "fuel_dens",
                                  [0., np.amax(fuel_dens[0][::1, ::1, iplane-1], axis=None)], ps.img_specs,
                                  ui=None, draw_arrow=False, fuel_green=True, color_burn = True)
            del fuel_dens
        return

def emissions(flags,qf,ps,plane=1):
    if flags.emissions == 2 or flags.emissions == 3:
        print("\t-pm emissions")
        emiss = read_fireca_field("pm_emissions-", qf.ntimes_ave, qf.time_ave, qf, 0, ps.output_folder)
        img.plot_2d_field(True, qf, plane, 'xy', emiss, "PM2.5 [g]", "pm_emissions",
                           [], ps.img_specs, ps.fuel_idx, ps.no_fuel_idx, flags)
        del emiss

    if flags.emissions == 1 or flags.emissions == 3:
        print("\t-co emissions")
        emiss = read_fireca_field("co_emissions-", qf.ntimes_ave, qf.time_ave, qf, 0, ps.output_folder)
        img.plot_2d_field(True, qf, plane, 'xy', emiss, "CO [g]", "co_emissions",
                           [], ps.img_specs, ps.fuel_idx, ps.no_fuel_idx, flags)
        del emiss

    if flags.emissions == 5:
        print("\t-water emissions")
        emiss = read_fireca_field("water_emissions-", qf.ntimes_ave, qf.time_ave, qf, 0, ps.output_folder)
        img.plot_2d_field(True, qf, plane, 'xy', emiss, "Water [g]", "water_emissions",
                           [], ps.img_specs, ps.fuel_idx, ps.no_fuel_idx, flags)
        del emiss
    return

def radiation(flags,qf,ps,plane=1):
    if flags.thermal_rad == 1:
        print("\t-radiation")
        conv = read_fireca_field("thermalradiation-", qf.ntimes_ave, qf.time_ave, qf, 0, ps.output_folder)
        img.plot_2d_field(True, qf, plane, 'xy', conv, "Convective heat to human [kW/m^2 skin]", "conv_heat",
                           [], ps.img_specs, ps.fuel_idx, ps.no_fuel_idx, flags)
        del conv
    return

def eng_to_atmos(flags,qf,ps,plane=1):
    if flags.en2atm == 1:
        print("\t-energy to atmosphere")
        en_to_atm = read_fireca_field("fire-energy_to_atmos-", qf.ntimes, qf.time, qf, 1, ps.output_folder)
        img.plot_2d_field(False, qf, plane, 'xy', en_to_atm, "Energy to atmosphere [kW/m^3]", "en_to_atm",
                           [], ps.img_specs, ps.fuel_idx, ps.no_fuel_idx, flags)
        del en_to_atm
    return

def qf_winds(flags,qf,ps,plane=1):
    if flags.qf_winds == 1:
        print("\t-fireca winds")
        print("    * u")
        windu_qf = read_fireca_field("windu", qf.ntimes, qf.time, qf, 1, ps.output_folder)
        img.plot_2d_field(False, qf, plane, 'xy', windu_qf, "U [m/s]", "u", [], ps.img_specs, [], [], flags)
        del windu_qf

        print("    * v")
        windv_qf = read_fireca_field("windv", qf.ntimes, qf.time, qf, 1, ps.output_folder)
        img.plot_2d_field(False, qf, plane, 'xy', windv_qf, "V [m/s]", "v", [], ps.img_specs, [], [], flags)
        del windv_qf

        print("    * w")
        windw_qf = read_fireca_field("windw", qf.ntimes, qf.time, qf, 1, ps.output_folder)
        img.plot_2d_field(False, qf, plane, 'xy', windw_qf, "W [m/s]", "w", [], ps.img_specs, [], [], flags)
        del windw_qf
    return

def qu_winds_inst(flags, qf, qu, ps, plane=2):
    if flags.qu_qwinds_inst == 1:
        vertplane = math.floor(qu.Ly*0.5 / qu.dy)
        print("\t-QU winds")
        print("    * u")
        windu_qu = read_fireca_field("qu_windu", qu.ntimes, qu.time, qu, 1, ps.output_folder)
        img.plot_2d_field(False, qu, plane, 'xy', windu_qu, "U_inst [m/s]", "u_qu", [], ps.img_specs, [], [], flags)
        img.plot_2dvert_field(False, qu, vertplane, 'xz', windu_qu, "U_inst [m/s]",
                          "u_qu_vert", [], ps.img_specs, [], flags)
        del windu_qu

        print("    * v")
        windv_qu = read_fireca_field("qu_windv", qu.ntimes, qu.time, qu, 1, ps.output_folder)
        img.plot_2d_field(False, qu, plane, 'xy', windv_qu, "V_inst [vm/s]", "v_qu", [], ps.img_specs, [], [], flags)
        img.plot_2dvert_field(False, qu, vertplane, 'xz', windv_qu, "V_inst [m/s]", "v_qu_vert",
                          [], ps.img_specs, [], flags)
        del windv_qu

        print("    * w")
        windw_qu = read_fireca_field("qu_windw", qu.ntimes, qu.time, qu, 1, ps.output_folder)
        img.plot_2d_field(False, qu, plane, 'xy', windw_qu, "W_inst [m/s]", "w_qu", [], ps.img_specs, [], [], flags)
        img.plot_2dvert_field(False, qu, vertplane, 'xz', windw_qu, "W_inst [m/s]", "w_qu_vert",
                          [], ps.img_specs, [], flags)
        img.plot_isosurface("qu_windw", qf, qu, qu.ntimes, windw_qu, ps.img_specs, flags, ps.output_folder)
        del windw_qu
    return

def qu_winds_avg(flags,qf,qu,ps,plane=2):
    if flags.qu_qwinds_ave == 1:
        print("\t-QU winds ave")
        print("    * u")
        windu_qu = read_fireca_field("qu_windu_ave", qu.ntimes_ave, qu.time_ave, qu, 1, ps.output_folder)
        img.plot_2d_field(True, qu, plane, 'xy', windu_qu, "U_ave [m/s]", "u_qu_ave", [], ps.img_specs, [], [], flags)
        del windu_qu

        print("    * v")
        windv_qu = read.read_fireca_field("qu_windv_ave", qu.ntimes_ave, qu.time_ave, qu, 1, ps.output_folder)
        img.plot_2d_field(True, qu, plane, 'xy', windv_qu, "V_ave [m/s]", "v_qu_ave", [], ps.img_specs, [], [], flags)
        del windv_qu

        print("    * w")
        windw_qu = read.read_fireca_field("qu_windw_ave", qu.ntimes_ave, qu.time_ave, qu, 1, ps.output_folder)
        img.plot_2d_field(True, qu, plane, 'xy', windw_qu, "W_ave [m/s]", "w_qu_ave", [], ps.img_specs, [], [], flags)
        img.plot_isosurface("qu_windw_ave", qf, qu, qu.ntimes, windw_qu, ps.img_specs, flags, ps.output_folder)
        del windw_qu
    return

def w_plumes(flags,qf,qu,ps,plane=2):
    if flags.qu_qwinds_inst == 1:
        print("\t-QU w-plumes")
        q = copy.deepcopy(qu)
        if qu.time[0] == 0:
            q.ntimes = qu.ntimes - 1
            q.times = qu.time[1:]

        wplume = read_fireca_field("qu_wplume", q.ntimes, q.times, q, 1, ps.output_folder)
        img.plot_2d_field(False, q, plane, 'xy', wplume, "W_plumes [m/s]", "wplume", [], ps.img_specs, [], [], flags)
    return

def reaction_rate(flags,qf,ps,plane=1):
    if flags.react_rate == 1:
        print("\t-reaction rate")
        react_rate = read_fireca_field("fire-reaction_rate-", qf.ntimes, qf.time, qf, 0, ps.output_folder)
        img.plot_2d_field(False, qf, plane, 'xy', react_rate,
                      "Reaction rate [kg/m3/s]", "react_rate", [], ps.img_specs, ps.fuel_idx, ps.no_fuel_idx, flags)
        del react_rate
    return

def fuel_moisture(flags, qf, ps):   
    if flags.moisture == 1:
        fuels_moist = read_fireca_field("fuels-moist-", qf.ntimes, qf.time, qf, 0, ps.output_folder)
        for iplane in ps.moist_planes:
            if iplane <= qf.nz:
                print("\t-moisture, plane: %d" % iplane)
                plot_2d_field(False, qf, iplane, fuels_moist,
                              "Fuel moisture [-]", "fuels_moist", [], ps.img_specs, 
                              ui=None, draw_arrow=False, fuel_green=True, color_burn = True)
        del fuels_moist
    return

def plot_2d_field(is_ave_time, q, plane, plotvar, ystr,
                  savestr, cblim, img_specs, 
                  ui=None, draw_arrow=False, fuel_green=False, color_burn = False):

    file_list = []
    
    if is_ave_time is True:
        ntimes = q.ntimes_ave
        times = q.time_ave
    else:
        ntimes = q.ntimes
        times = q.time

    if not cblim:
        myvmin = 1e8
        myvmax = -1e8
        for i in range(0, ntimes):
            if plane:
                currval = plotvar[i][::1, ::1, plane - 1]
            else:
                currval = plotvar[i]

            myvmin = min(myvmin, np.amin(currval))
            myvmax = max(myvmax, np.amax(currval))

        myvmin = math.floor(myvmin)
        myvmax = math.ceil(myvmax)
        if myvmin == myvmax:
            if myvmin == 0:
                myvmin = -1
                myvmax = +1
            else:
                myvmin *= 0.5
                myvmax *= 2.
    else:
        myvmin = cblim[0]
        myvmax = cblim[1]

    for i in range(0, ntimes):
        #print("     * time %d/%d" % (i + 1, ntimes))

        if plane:
            if plane<0:
                #ZC take z indexes from UI for midstory and canopy figures
                if plane== ui.ms_bottom:
                    pmax = ui.ms_top
                else:
                    pmax = plotvar[0].shape[-1]
                refval = np.sum( plotvar[0][::1, ::1, abs(plane) - 1:pmax:1], axis=2 ) +1.e-6
                currval = (np.sum( plotvar[i][::1, ::1, abs(plane) - 1:pmax:1], axis=2 ))/refval
                currval = np.where(refval>1.e-6, currval, -1.)
            else:
                currval = plotvar[i][::1, ::1, plane - 1]
            plane_str = '_Plane_%d' % abs(plane)
        else:
            currval = plotvar[i]
            plane_str = ''

        currval = currval.squeeze()

        fig = pylab.figure(figsize=(img_specs.figure_size[0], img_specs.figure_size[1]))
        ax = fig.add_subplot(111)
        
        if fuel_green:
            pylab.imshow(currval, cmap='YlGn', interpolation='none', origin='lower',
                         extent=q.horizontal_extent, vmin=myvmin, vmax=myvmax)
        else: 
            pylab.imshow(currval, cmap='jet', interpolation='none', origin='lower',
                         extent=q.horizontal_extent, vmin=myvmin, vmax=myvmax)
        cbar = pylab.colorbar()
        cbar.set_label(ystr, size=img_specs.axis_font["size"], fontname=img_specs.axis_font["fontname"])
        cbar.ax.tick_params(labelsize=img_specs.colorbar_font["size"])
        pylab.xlabel('X [m]', **img_specs.axis_font)
        pylab.ylabel('Y [m]', **img_specs.axis_font)
        pylab.title('Time = %s s' % times[i], **img_specs.title_font)
        set_ticks_font(img_specs.axis_font, ax)
        
        #ZC draw wind direction arrow with speed label for gif
        if draw_arrow==True:
            
            #length of domain in meters            
            x_len, y_len = q.horizontal_extent[1], q.horizontal_extent[3]
            
            #calculated length and width to determine where to place arrow and txt box
            if y_len>x_len:
                length = x_len   
            else:
                length = y_len 
            
            txt = pylab.text(length*0.95, length*0.05, str(int(q.avg_wind_speeds[i]))+"\n[m/s]",
                        fontsize= img_specs.axis_font['size'], backgroundcolor="white", c="black", ma="center", ha="right", 
                        va="bottom")
            
            pylab.draw()
            patch = txt.get_bbox_patch()
            box  = patch.get_extents()
            tcbox = ax.transData.inverted().transform(box)
            
            start_x = np.average(tcbox[:,0])
            arrow_length = tcbox[1,0]-start_x
            width = arrow_length / 4
            start_y = tcbox[1,1] + arrow_length*1.25
            
            #Draw arrow    
            end_x, end_y = pol2cart(arrow_length, q.avg_wind_directions[i])
            pylab.arrow(start_x, start_y, end_x, end_y, fc="white", ec="black", 
                        shape='full', width=width, head_width=width*3, 
                        head_length=width*3, length_includes_head=True)
            
        if color_burn:
            ##Color burned and burning cells
            fire_locations = np.where(plotvar[0][:,:,0]!=plotvar[i][:,:,0])
            plot_color(fire_locations, currval, 'red', myextent=q.horizontal_extent)
            
            burned_locations = np.where((plotvar[0][:,:,0]!=plotvar[i][:,:,0]) & (plotvar[i][:,:,0]< 0.01))
            plot_color(burned_locations, currval, 'black', myextent=q.horizontal_extent)
        
        fname = '%s_Time_%d_s%s.png' % (savestr, times[i], plane_str)
        fname = os.path.join(img_specs.save_dir, fname)
        pylab.savefig(fname)
        pylab.close()
        if img_specs.gen_gif == 1:
            file_list.append(fname)

    if img_specs.gen_gif == 1:
        fname_gif = '%s%s.gif' % (savestr, plane_str)
        fname_gif = os.path.join(img_specs.gif_dir, fname_gif)
        img.make_gif(fname_gif, file_list)

def set_ticks_font(axis_font: dict, ax):
    # https://stackoverflow.com/questions/3899980/how-to-change-the-font-size-on-a-matplotlib-plot
    for label in (ax.get_xticklabels() + ax.get_yticklabels()):
        label.set_fontname(axis_font["fontname"])
        label.set_fontsize(axis_font["size"])

    try:
        for label in ax.get_zticklabels():
            label.set_fontname(axis_font["fontname"])
            label.set_fontsize(axis_font["size"])
    except:
        return

def plot_color(locations, currval, color, myextent):
    plot_bool = np.zeros_like(currval)
    plot_bool[locations] = 1
    #Build cmap
    cmap = mpl.colors.LinearSegmentedColormap.from_list('my_cmap',['white', color],256)
    cmap._init() # create the _lut array, with rgba values
    # create your alpha array and fill the colormap with them.
    # here it is progressive, but you can create whathever you want
    alphas = np.linspace(0, 0.9, cmap.N+3)
    cmap._lut[:,-1] = alphas
    
    pylab.imshow(plot_bool,cmap=cmap, interpolation='none', origin='lower',
                 extent=myextent, vmin=0, vmax=1)

def pol2cart(rho, phi):
                """
                ZC. This functions converts polar cordinates to cartisian coordinates
                
                Inputs:
                    rho: magnitude (float/int) 
                    phi: direction (float/int)

                Returns
                    x and y index
                """
                x = -rho * np.sin(np.radians(phi))
                y = -rho * np.cos(np.radians(phi))
                return(x, y)

if __name__ == '__main__':
    # Parameters:
    # 1) project folder. Defaults to working directory. Use "default path" if defining further arguments
    # 2) planes to visualize fuel mositure ()
    # 3) planes to visualize fule density ()
    
    if len(sys.argv) == 1:
        prj_folder_in = os.getcwd()
    else:
        prj_folder_in = sys.argv[1]


    main(prj_folder_in)