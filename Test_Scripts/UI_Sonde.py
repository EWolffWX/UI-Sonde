import tkinter as tk
from tkinter import ttk
from tkinter import PhotoImage
#from tkinter.filedialog import askopenfilename
from tkinter import filedialog
from tkinter import *

import matplotlib.pyplot as plt
import matplotlib.axes as axes
import pandas as pd
import math
import argparse
import os

import metpy.calc as mpcalc
from metpy.cbook import get_test_data
from metpy.plots import add_metpy_logo, SkewT, Hodograph
from metpy.units import units
from metpy.plots import USCOUNTIES
import metpy.interpolate as minterp

import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.mpl.geoaxes as geo

#import pyart
#import fsspec
#from geopy.distance import geodesic

from PIL import Image


def plot_sonde(inputlat, inputlon, inputfile, inputprof, inputdpi, inputlogo):

    station_lat = float(inputlat)
    station_lon = float(inputlon)
    user_filepath = str(inputfile)
    prof_choice = str(inputprof)
    dpi_choice = int(inputdpi)
    user_logo = str(inputlogo)

    
    # Define functions
    def find_temp(pressure_value, offset):
        a = 80.18514962
        b = -240.03647482 + offset

        y_fit = (a * np.log10(pressure_value.magnitude)) + b 

        return y_fit

    def hgt_fit(pressure_values, height_values, pressure_point):
        x = p.magnitude
        y = hgt.magnitude
        xlog = np.log10(x)
        curve = np.polyfit(xlog, y, 1)
        y_fit = ((curve[0] * np.log10(pressure_point.magnitude)) + curve[1])*units.meters
        return y_fit

    # Courtesy MetPy
    def effective_layer(p, t, td, h, height_layer=False):
        """A function that determines the effective inflow layer for a convective sounding.

        Uses the default values of Thompason et al. (2004) for CAPE (100 J/kg) and CIN (-250 J/kg).

        Input:
          - p: sounding pressure with units
          - T: sounding temperature with units
          - Td: sounding dewpoint temperature with units
          - h: sounding heights with units

        Returns:
          - pbot/hbot, ptop/htop: pressure/height of the bottom level,
                                  pressure/height of the top level
        """
        from metpy.calc import cape_cin, parcel_profile
        from metpy.units import units

        pbot = None

        for i in range(p.shape[0]):
            prof = parcel_profile(p[i:], t[i], td[i])
            sbcape, sbcin = cape_cin(p[i:], t[i:], td[i:], prof)
            if sbcape >= 100 * units('J/kg') and sbcin > -250 * units('J/kg'):
                pbot = p[i]
                hbot = h[i]
                bot_idx = i
                break
        if not pbot:
            return None, None

        for i in range(bot_idx + 1, p.shape[0]):
            prof = parcel_profile(p[i:], t[i], td[i])
            sbcape, sbcin = cape_cin(p[i:], t[i:], td[i:], prof)
            if sbcape < 100 * units('J/kg') or sbcin < -250 * units('J/kg'):
                ptop = p[i]
                htop = h[i]
                break

        if height_layer:
            return hbot, htop
        else:
            return pbot, ptop

    def wind_damage(MLCAPE, MLCIN, LR03, Mean_Wind):
        """Function to compute the wind damage parameter for a convective sounding.

        Modified from the SHARPpy package.

        Input:
          - MLCAPE: mixed-layer CAPE (J/kg)
          - MLCIN: mixed-layer CIN (J/kg)
          - LR03: 0-3 km lapse rate (must be unitless)
          - Mean_Wind: 1-3.5 km mean wind (needs to be in m/s)

        Returns:
          - wndg: wind damage parameter (unitless)
        """

        if LR03 < 7:
            LR03 = 0

        if MLCIN < -50*units('J/kg'):
            MLCIN = -50*units('J/kg')

        wndg = ((MLCAPE / 2000*units('J/kg')) * (LR03 / 9) * (Mean_Wind / 15*units('m/s')) * 
                ((50*units('J/kg') + MLCIN)/40*units('J/kg')))
        return wndg

    def ship_calc(mu_cape, mu_p, mu_T, mu_Td, p, T, shear06):
        '''
        Calculate the Sig Hail Parameter (SHIP)

        Modified from the SHARPpy package.
        
        Ryan Jewell (SPC) helped in correcting this equation as the SPC
        sounding help page version did not have the correct information
        of how SHIP was calculated.

        The significant hail parameter (SHIP; SPC 2014) is
        an index developed in-house at the SPC. (Johnson and Sugden 2014)

        Parameters
        ----------
        
        mucape: pint quantity
            Most unstable CAPE (J/kg)
        mu_p: pint quantity
            Most unstable parcel pressure (hPA)
        mu_T: pint quantity
            Most unstable parcel temperature (C)
        mu_Td: pint quantity
            Most unstable parcel dew point (C)
        p: pint quantity
            Pressure profile (hPA)
        T: pint quantity
            Temperature profile (C)
        shear06: pint quantity
            0-6 km shear (m/s)

        Returns
        -------
        ship : number
            significant hail parameter (unitless)

        '''

        # Compute LR75
        T500 = minterp.interpolate_1d(500*units.hPa, p, T)
        T700 = minterp.interpolate_1d(700*units.hPa, p, T)
        H500 = (minterp.interpolate_1d(500*units.hPa, p, hgt)).to(units.km)
        H700 = (minterp.interpolate_1d(700*units.hPa, p, hgt)).to(units.km)
        LR75 = ((T500 - T700)/(H700-H500)).magnitude

        # Compute Freezing Level
        frz_lvl = minterp.interpolate_1d(0, T.magnitude, hgt.magnitude) * units.m

        # Compute MU parcel mixing ratio
        mu_rh = mpcalc.relative_humidity_from_dewpoint(mu_T, mu_Td)
        mu_mr = (mpcalc.mixing_ratio_from_relative_humidity(mu_p, mu_T, mu_rh))*1000


        shear06 = shear06.magnitude

        if shear06 > 27:
            shear06 = 27.
        elif shear06 < 7:
            shear06 = 7.

        if mu_mr > 13.6:
            mu_mr = 13.6
        elif mu_mr < 11.:
            mu_mr = 11.

        if T500.magnitude > -5.5:
            T500.magnitude = -5.5

        ship = -1. * (mu_cape.magnitude * mu_mr * LR75 * T500.magnitude * shear06) / 42000000.

        if mu_cape.magnitude < 1300:
            ship = ship*(mu_cape.magnitude/1300.)

        if LR75 < 5.8:
            ship = ship*(LR75/5.8)

        if frz_lvl.magnitude < 2400:
            ship = ship * (frz_lvl.magnitude/2400.)

        return ship[0]
    
    
    
    # Parse input data

    launch_lat = station_lat
    launch_lon = station_lon

    filepath = user_filepath

    
    # Open text file
    df_header = pd.read_table(filepath, delimiter=',',  
                     skiprows=0,
                     names=['LEV', 'HGHT', 'T', 'TD', 'WDIR', 'WSPD']
                     )
    head_info = df_header['LEV'][1].split()
    crew_name = head_info[0]
    if crew_name == 'XXX':
        crew_name = '---'
    year = int(head_info[1][:2])
    month = int(head_info[1][2:4])
    day = int(head_info[1][4:6])
    hour = int(head_info[1][7:9])
    minute = int(head_info[1][9:11])

    # Format launch time
    if year <= 9:
        year = '200' + str(year)
    else:
        year = '20'+ str(year)
    if month <= 9:
        month = '0' + str(month)
    else:
        month = str(month)
    if day <= 9:
        day = '0' + str(day)
    else:
        day = str(day)
    if hour <= 9:
        hour = '0' + str(hour)
    else:
        hour = str(hour)
    if minute <= 9:
        minute = '0' + str(minute)
    else:
        minute = str(minute)

    df = pd.read_table(filepath, delimiter=',',  
                     skiprows=6,
                     names=['LEV', 'HGHT', 'T', 'TD', 'WDIR', 'WSPD']
                     )
    # Drop any rows with all NaN values for T, Td, winds
    df = df.dropna(subset=('T', 'TD', 'WDIR', 'WSPD'), how='all'
                   ).reset_index(drop=True)
    df = df[df['T'] != -9999.0].reset_index(drop=True)
    df = df[df['WSPD'] != -9999.0].reset_index(drop=True)

    new_p = []
    for i in np.arange(0, len(df['LEV']), 1):
        temp_lev = float(df['LEV'][i])
        new_p.append(temp_lev)

    # Remove erroneous surface?
    remove_thresh = 0 #Number of initial values to remove from sounding

    # Add Metpy units
    p = new_p[remove_thresh:] * units.hPa
    T = df['T'].values[remove_thresh:] * units.degC
    Td = df['TD'].values[remove_thresh:] * units.degC
    wind_speed = df['WSPD'].values[remove_thresh:] * units.knots
    wind_dir = df['WDIR'].values[remove_thresh:] * units.degrees
    u, v = mpcalc.wind_components(wind_speed, wind_dir)
    hgt = df['HGHT'].values[remove_thresh:] * units.meter

    # Calculate parcel profiles and related params
    # Sfc
    lcl_pressure, lcl_temperature = mpcalc.lcl(p[0], T[0], Td[0])
    el_pressure, _ = mpcalc.el(p, T, Td)
    prof = mpcalc.parcel_profile(p, T[0], Td[0]).to('degC')
    lfc_pressure, _ = mpcalc.lfc(p, T, Td, prof, which='most_cape')
    #lcl_h = hgt_fit(p, hgt, lcl_pressure)
    lcl_h = minterp.interpolate_1d(lcl_pressure, p, hgt)
    el_h = minterp.interpolate_1d(el_pressure, p, hgt)
    lfc_h = minterp.interpolate_1d(lfc_pressure, p, hgt)
    sfc_cape, sfc_cin = mpcalc.surface_based_cape_cin(p, T, Td)
    # ML
    ml_p, ml_T, ml_Td = mpcalc.mixed_parcel(p, T, Td, depth = 100 * units.hPa)
    ml_prof = mpcalc.parcel_profile(p, ml_T, ml_Td).to('degC')
    ml_lcl_pressure, ml_lcl_temperature = mpcalc.lcl(ml_p, ml_T, ml_Td)
    ml_el_pressure, _ = mpcalc.el(p, T, Td, parcel_temperature_profile = ml_prof)
    ml_lfc_pressure, _ = mpcalc.lfc(p, T, Td, ml_prof, which='most_cape')
    ml_lcl_h = minterp.interpolate_1d(ml_lcl_pressure, p, hgt)
    ml_el_h = minterp.interpolate_1d(ml_el_pressure, p, hgt)
    ml_lfc_h = minterp.interpolate_1d(ml_lfc_pressure, p, hgt)
    ml_cape, ml_cin = mpcalc.mixed_layer_cape_cin(p, T, Td)
    # MU
    mu_p, mu_T, mu_Td, mu_ind = mpcalc.most_unstable_parcel(p, T, Td)
    mu_prof = mpcalc.parcel_profile(p[mu_ind:], mu_T, mu_Td).to('degC')
    mu_lcl_pressure, mu_lcl_temperature = mpcalc.lcl(mu_p, mu_T, mu_Td)
    mu_el_pressure, _ = mpcalc.el(p[mu_ind:], T[mu_ind:], Td[mu_ind:], parcel_temperature_profile = mu_prof)
    mu_lfc_pressure, _ = mpcalc.lfc(p[mu_ind:], T[mu_ind:], Td[mu_ind:], mu_prof, which='most_cape')
    mu_lcl_h = minterp.interpolate_1d(mu_lcl_pressure, p, hgt)
    mu_el_h = minterp.interpolate_1d(mu_el_pressure, p, hgt)
    mu_lfc_h = minterp.interpolate_1d(mu_lfc_pressure, p, hgt)
    mu_cape, mu_cin = mpcalc.most_unstable_cape_cin(p, T, Td)

    # Calculate additional parameters
    pw = mpcalc.precipitable_water(p, Td).to(units.inches)
    if ((p[-1:][0] < 500*units.hPa) & (sfc_cape >= 10*units('J/kg'))):
        dcape, _, _ = mpcalc.downdraft_cape(p, T, Td)
    else:
        dcape = np.nan * units('J/kg')

    if hgt[-1:][0] < (6000*units.meters+hgt[0]):
        rm = [0 * units.knots, 0 * units.knots]
        lm = [0 * units.knots, 0 * units.knots]
        mean_storm = [0 * units.knots, 0 * units.knots]
    else:
        rm, lm, mean_storm = mpcalc.bunkers_storm_motion(p, u, v, hgt)
    # Shear and SRH
    if hgt[-1:][0] < (1000*units.meters+hgt[0]):
        shear01 = np.nan * units.knots
        srh01 = np.nan * units('m/s^2')
    else:
        u01, v01 = mpcalc.bulk_shear(p, u, v, height = hgt, depth = 1000*units.meters, bottom = hgt[0])
        shear01 = np.sqrt((u01**2) + (v01**2))
        srh01,_,_ = mpcalc.storm_relative_helicity(hgt, u, v, depth=1000*units.meters, storm_u=rm[0], storm_v=rm[1])

    if hgt[-1:][0] < (3000*units.meters+hgt[0]):
        shear03 = np.nan * units.knots
        srh03 = np.nan * units('m/s^2')
    else:
        u03, v03 = mpcalc.bulk_shear(p, u, v, height = hgt, depth = 3000*units.meters, bottom = hgt[0])
        shear03 = np.sqrt((u03**2) + (v03**2))
        srh03,_,_ = mpcalc.storm_relative_helicity(hgt, u, v, depth=3000*units.meters, storm_u=rm[0], storm_v=rm[1])

    if hgt[-1:][0] < (6000*units.meters+hgt[0]):
        shear06 = np.nan * units.knots
    else:
        u06, v06 = mpcalc.bulk_shear(p, u, v, height = hgt, depth = 6000*units.meters, bottom = hgt[0])
        shear06 = np.sqrt((u06**2) + (v06**2))

    if hgt[-1:][0] < (10000*units.meters+hgt[0]):
        shear510 = np.nan * units.knots
    else:
        u510, v510 = mpcalc.bulk_shear(p, u, v, height = hgt, depth = 5000*units.meters, bottom = 5000*units.meters)
        shear510 = np.sqrt((u510**2) + (v510**2))
    # Effective layer stuff
    pbot, ptop = effective_layer(p, T, Td, hgt, height_layer=False)
    hbot, htop = effective_layer(p, T, Td, hgt, height_layer=True)
    if hbot is not None:
        ueff, veff = mpcalc.bulk_shear(p, u, v, height = hgt, depth = (htop-hbot), bottom = hbot)
        sheareff = np.sqrt((ueff**2) + (veff**2))
        srheff,_,_ = mpcalc.storm_relative_helicity(hgt, u, v, depth= (htop-hbot), bottom = hbot, 
                                                    storm_u=rm[0], storm_v=rm[1])
    else:
        sheareff = np.nan * units.knots
        srheff = np.nan * units('m/s^2')
    # Composites
    if (not math.isnan(srheff.magnitude)) & (not math.isnan(sheareff.magnitude)):
        supercell = mpcalc.supercell_composite(mu_cape, srheff, sheareff)
    else:
        supercell = np.nan * units.dimensionless
    if (not math.isnan(srh01.magnitude)) & (not math.isnan(shear06.magnitude)):  
        stp = mpcalc.significant_tornado(sfc_cape, lcl_h, srh01, shear06.to(units('m/s')))
    else:
        stp = np.nan * units.dimensionless
    # Extra params for hazard type
    ind_1km = np.where(hgt>=1000*units.meters)[0][0]
    LR01 = -1 * ((T[ind_1km]-T[0]).magnitude)
    if hgt[-1:][0] < (8000*units.meters+hgt[0]):
        shear08 = np.nan * units.knots
    else:
        u08, v08 = mpcalc.bulk_shear(p, u, v, height = hgt, depth = 8000*units.meters, bottom = hgt[0])
        shear08 = np.sqrt((u08**2) + (v08**2))
    #
    if hgt[-1:][0] < (3500*units.meters+hgt[0]):
        wndg = np.nan * units.dimensionless
    else:
        ind_3km = np.where(hgt>=3000*units.meters)[0][0]
        LR03 = -1/3 * ((T[ind_3km]-T[0]).magnitude)
        ind_1km = np.where(hgt>=1000*units.meters)[0][0]
        ind_3_5km = np.where(hgt>=3500*units.meters)[0][0]

        Mean_Wind = np.mean(mpcalc.wind_speed(u[ind_1km:ind_3_5km], v[ind_1km:ind_3_5km]).to(units('m/s')))
        wndg = (wind_damage(ml_cape, ml_cin, LR03, Mean_Wind).magnitude)*units.dimensionless
        
    if (p[-1:][0] > 700*units.hPa) | (hgt[-1:][0] < 6000*units.meters):
        ship = np.nan * units.dimensionless
    else:
        ship = ship_calc(mu_cape, mu_p, mu_T, mu_Td, p, T, shear06)
    srf_RH = mpcalc.relative_humidity_from_dewpoint(T[0], Td[0])
    srf_HI = mpcalc.heat_index(T[0], srf_RH)[0].to(units.degF)
    #print(srf_HI.magnitude)




    #### Plotting time! ####

    # Make figure
    fig = plt.figure(figsize=(9, 9), facecolor='k', constrained_layout=True)
    #fig, ax = plt.subplots(figsize=(9, 9), constrained_layout=True) 

    # Plot Skew-T
    skew = SkewT(fig, rotation=45)
    skew.plot(p, T, 'r', linewidth=2)
    skew.plot(p, Td, 'limegreen', linewidth=2)
    barb_loc = int(len(df)/30)
    skew.plot_barbs(p[::barb_loc], u[::barb_loc], v[::barb_loc], barbcolor='w')
    skew.ax.set_ylim(1000, 90)
    skew.ax.set_xlim(-50, 50)

    skew.ax.set_xlabel(f'Temperature ({T.units:~P})', fontsize=12, c='w')
    skew.ax.set_ylabel(f'Pressure ({p.units:~P})', fontsize=12, c='w')
    skew.ax.tick_params(labelsize=12, labelcolor='w', color='w')
    skew.ax.set_facecolor('k')

    # Plot parcel profile and shade CAPE/CIN
    # Add lines for lcl, lfc, el
    if prof_choice == 'Surface':
        skew.plot(p, prof, 'w', linewidth=2, linestyle='--')
        skew.shade_cin(p, T, prof, Td)
        skew.shade_cape(p, T, prof)
        skew.plot(np.array([lcl_pressure.magnitude, lcl_pressure.magnitude])*units.hPa, 
                  np.array([(find_temp(lcl_pressure, 30)), (find_temp(lcl_pressure, 30)+5)])*units.degC,
                  'limegreen', linewidth=2)
        skew.ax.text((find_temp(lcl_pressure, 30)+5.5), lcl_pressure.magnitude, 
                     'LCL', c='limegreen', va='center', ha='left', fontsize=13, fontweight='bold')
        skew.plot(np.array([el_pressure.magnitude, el_pressure.magnitude])*units.hPa, 
                  np.array([(find_temp(el_pressure, 30)), (find_temp(el_pressure, 30)+5)])*units.degC,
                  'magenta', linewidth=2)
        skew.ax.text((find_temp(el_pressure, 30)+5.5), el_pressure.magnitude, 
                     'EL', c='magenta', va='center', ha='left', fontsize=13, fontweight='bold')
        skew.plot(np.array([lfc_pressure.magnitude, lfc_pressure.magnitude])*units.hPa, 
                  np.array([(find_temp(lfc_pressure, 30)), (find_temp(lfc_pressure, 30)+5)])*units.degC,
                  'yellow', linewidth=2)
        skew.ax.text((find_temp(lfc_pressure, 30)+5.5), lfc_pressure.magnitude, 
                     'LFC', c='yellow', va='center', ha='left', fontsize=13, fontweight='bold')
    elif prof_choice == 'Mixed Layer':
        skew.plot(p, ml_prof, 'w', linewidth=2, linestyle='--')
        skew.shade_cin(p, T, ml_prof, Td)
        skew.shade_cape(p, T, ml_prof)
        skew.plot(np.array([ml_lcl_pressure.magnitude, ml_lcl_pressure.magnitude])*units.hPa, 
                  np.array([(find_temp(ml_lcl_pressure, 30)), (find_temp(ml_lcl_pressure, 30)+5)])*units.degC,
                  'limegreen', linewidth=2)
        skew.ax.text((find_temp(ml_lcl_pressure, 30)+5.5), ml_lcl_pressure.magnitude, 
                     'LCL', c='limegreen', va='center', ha='left', fontsize=13, fontweight='bold')
        skew.plot(np.array([ml_el_pressure.magnitude, ml_el_pressure.magnitude])*units.hPa, 
                  np.array([(find_temp(ml_el_pressure, 30)), (find_temp(ml_el_pressure, 30)+5)])*units.degC,
                  'magenta', linewidth=2)
        skew.ax.text((find_temp(ml_el_pressure, 30)+5.5), ml_el_pressure.magnitude, 
                     'EL', c='magenta', va='center', ha='left', fontsize=13, fontweight='bold')
        skew.plot(np.array([ml_lfc_pressure.magnitude, ml_lfc_pressure.magnitude])*units.hPa, 
                  np.array([(find_temp(ml_lfc_pressure, 30)), (find_temp(ml_lfc_pressure, 30)+5)])*units.degC,
                  'yellow', linewidth=2)
        skew.ax.text((find_temp(ml_lfc_pressure, 30)+5.5), ml_lfc_pressure.magnitude, 
                     'LFC', c='yellow', va='center', ha='left', fontsize=13, fontweight='bold')
    elif prof_choice == 'Most Unstable':
        skew.plot(p[mu_ind:], mu_prof, 'w', linewidth=2, linestyle='--')
        skew.shade_cin(p[mu_ind:], T[mu_ind:], mu_prof, Td[mu_ind:])
        skew.shade_cape(p[mu_ind:], T[mu_ind:], mu_prof)
        skew.plot(np.array([mu_lcl_pressure.magnitude, mu_lcl_pressure.magnitude])*units.hPa, 
                  np.array([(find_temp(mu_lcl_pressure, 30)), (find_temp(mu_lcl_pressure, 30)+5)])*units.degC,
                  'limegreen', linewidth=2)
        skew.ax.text((find_temp(mu_lcl_pressure, 30)+5.5), mu_lcl_pressure.magnitude, 
                     'LCL', c='limegreen', va='center', ha='left', fontsize=13, fontweight='bold')
        skew.plot(np.array([mu_el_pressure.magnitude, mu_el_pressure.magnitude])*units.hPa, 
                  np.array([(find_temp(mu_el_pressure, 30)), (find_temp(mu_el_pressure, 30)+5)])*units.degC,
                  'magenta', linewidth=2)
        skew.ax.text((find_temp(mu_el_pressure, 30)+5.5), mu_el_pressure.magnitude, 
                     'EL', c='magenta', va='center', ha='left', fontsize=13, fontweight='bold')
        skew.plot(np.array([mu_lfc_pressure.magnitude, mu_lfc_pressure.magnitude])*units.hPa, 
                  np.array([(find_temp(mu_lfc_pressure, 30)), (find_temp(mu_lfc_pressure, 30)+5)])*units.degC,
                  'yellow', linewidth=2)
        skew.ax.text((find_temp(mu_lfc_pressure, 30)+5.5), mu_lfc_pressure.magnitude, 
                     'LFC', c='yellow', va='center', ha='left', fontsize=13, fontweight='bold')

    # Add lines for effective
    if pbot is not None:
        skew.plot(np.array([pbot.magnitude, ptop.magnitude])*units.hPa,
                  np.array([(find_temp(pbot, -20)), (find_temp(ptop, -20))])*units.degC,
                  'cyan', linewidth=2)
        skew.plot(np.array([pbot.magnitude, pbot.magnitude])*units.hPa,
                  np.array([(find_temp(pbot, -22)), (find_temp(pbot, -18))])*units.degC,
                  'cyan', linewidth=2)
        skew.plot(np.array([ptop.magnitude, ptop.magnitude])*units.hPa,
                  np.array([(find_temp(ptop, -22)), (find_temp(ptop, -18))])*units.degC,
                  'cyan', linewidth=2)
        skew.ax.text((find_temp(ptop, -23)), ptop.magnitude, f'{int(htop.magnitude)} m', 
                     va='center', ha='right', fontsize=10, c='cyan')
        skew.ax.text((find_temp(pbot, -23)), pbot.magnitude, f'{int(hbot.magnitude)} m', 
                     va='center', ha='right', fontsize=10, c='cyan')


    # Add the relevant special lines
    skew.plot_dry_adiabats(linewidth=1)
    skew.plot_moist_adiabats(linewidth=1)
    skew.plot_mixing_lines()

    # Show the plot
    skew.ax.set_title(f'{crew_name} (Observed)', fontsize=18, loc='left', c='w', fontweight='bold')
    skew.ax.set_title(f'{year}-{month}-{day} {hour}:{minute} UTC', loc='right', fontsize=15, 
                       c='w', fontweight='bold')
    for spine in skew.ax.spines.values():
            spine.set_edgecolor('w')
    #ax.set_facecolor('k')


    # Plot hodograph
    cutoff_locs = [np.where(hgt.magnitude<3000)[0][-1], np.where(hgt.magnitude<6000)[0][-1],
                  np.where(hgt.magnitude<9000)[0][-1]]
    ax_hod = plt.axes([1.15, 0.40, 0.45, 0.45])
    h = Hodograph(ax_hod, component_range=70.)
    h.add_grid(increment=20)
    h.plot(u[0:cutoff_locs[0]], v[0:cutoff_locs[0]], c='r')  # Plot for 0-3 km
    if hgt[-1:][0] > 3000*units.meters:
        h.plot(u[cutoff_locs[0]:cutoff_locs[1]], v[cutoff_locs[0]:cutoff_locs[1]], 
               c='limegreen')  # Plot for 3-6 km
    if hgt[-1:][0] > 6000*units.meters:
        h.plot(u[cutoff_locs[1]-1:cutoff_locs[2]], v[cutoff_locs[1]-1:cutoff_locs[2]], 
               c='yellow')  # Plot for 6-9 km
    if hgt[-1:][0] > 9000*units.meters:
        h.plot(u[cutoff_locs[2]-1:], v[cutoff_locs[2]-1:], 
               c='cyan')  # Plot for 6-9 km
    h.ax.scatter(rm[0], rm[1], s=100, c='w', marker=9)
    h.ax.scatter(lm[0], lm[1], s=100, c='w', marker=8)

    if hbot is not None:
        hbot_ind = np.where(hgt.magnitude<=hbot.magnitude)[0][-1]
        htop_ind = np.where(hgt.magnitude<=htop.magnitude)[0][-1]
        h.ax.plot(np.array([rm[0].magnitude, u[hbot_ind].magnitude])*units.knots, 
                  np.array([rm[1].magnitude, v[hbot_ind].magnitude])*units.knots, c='cyan', linewidth=1)
        h.ax.plot(np.array([rm[0].magnitude, u[htop_ind].magnitude])*units.knots, 
                  np.array([rm[1].magnitude, v[htop_ind].magnitude])*units.knots, c='cyan', linewidth=1)

    ax_hod.set_xlabel('')
    ax_hod.set_ylabel('')
    ax_hod.tick_params(labelcolor='w', color='w')
    ax_hod.set_aspect(1)
    ax_hod.set_facecolor('k')
    for spine in ax_hod.spines.values():
            spine.set_edgecolor('w')
    # Set second title
    ax_hod.set_title('Plotted with UI-Sonde v0.4.0 + MetPy', loc='right', c='w')


    # Plot map
    ax_map = plt.axes([1.2, 0.06, 0.35, 0.35], projection=ccrs.PlateCarree())
    ax_map.add_feature(cfeature.STATES, linewidth=1.4, edgecolor='w')
    ax_map.add_feature(cfeature.LAND, facecolor='k')
    ax_map.add_feature(cfeature.OCEAN, facecolor='k')
    #ax_map.add_feature(USCOUNTIES, linewidth=0.5, edgecolor='gray')
    ax_map.set_extent([launch_lon-8, launch_lon+8, launch_lat-5, launch_lat+5])
    ax_map.scatter(launch_lon, launch_lat, color='r', s=300, marker='*', zorder=10)
    ax_map.set_aspect(1)
    for spine in ax_map.spines.values():
            spine.set_edgecolor('w')
   

    # Plot wind height
    ax_wind = plt.axes([0.98, 0.15, 0.1, 0.7])
    ax_wind.set_xticks([])  # Add x-axis ticks
    ax_wind.set_yticks([])  # Remove y-axis ticks
    #ax_wind.tick_params(labelsize=12, labelcolor='w', color='w')
    ax_wind.set_facecolor('k')
    for spine in ax_wind.spines.values():
            spine.set_edgecolor('w')
    ax_wind.set_xlim([0,140])
    ax_wind.set_ylim([1000,90])
    ax_wind.set_yscale('log')
    for w in np.arange(0, len(wind_speed.magnitude[0:cutoff_locs[0]])):
        wind_array = wind_speed.magnitude[0:cutoff_locs[0]]
        p_array = p.magnitude[0:cutoff_locs[0]]
        ax_wind.plot([0, wind_array[w]], [p_array[w], p_array[w]], c='r')
    for w in np.arange(0, len(wind_speed.magnitude[cutoff_locs[0]:cutoff_locs[1]])):
        wind_array = wind_speed.magnitude[cutoff_locs[0]:cutoff_locs[1]]
        p_array = p.magnitude[cutoff_locs[0]:cutoff_locs[1]]
        ax_wind.plot([0, wind_array[w]], [p_array[w], p_array[w]], c='limegreen')
    for w in np.arange(0, len(wind_speed.magnitude[cutoff_locs[1]:cutoff_locs[2]])):
        wind_array = wind_speed.magnitude[cutoff_locs[1]:cutoff_locs[2]]
        p_array = p.magnitude[cutoff_locs[1]:cutoff_locs[2]]
        ax_wind.plot([0, wind_array[w]], [p_array[w], p_array[w]], c='yellow')
    for w in np.arange(0, len(wind_speed.magnitude[cutoff_locs[2]:])):
        wind_array = wind_speed.magnitude[cutoff_locs[2]:]
        p_array = p.magnitude[cutoff_locs[2]:]
        ax_wind.plot([0, wind_array[w]], [p_array[w], p_array[w]], c='cyan')
    # Add ticks manually
    ax_wind.plot([40,40], [1000, 90], c='peru', linestyle='--')
    ax_wind.text(40, 1050, '40', c='w', ha='center', va='center')
    ax_wind.plot([80,80], [1000, 90], c='peru', linestyle='--')
    ax_wind.text(80, 1050, '80', c='w', ha='center', va='center')
    ax_wind.plot([120,120], [1000, 90], c='peru', linestyle='--')
    ax_wind.text(120, 1050, '120', c='w', ha='center', va='center')
    # Add title
    ax_wind.text(9, 102, 'Wind Speed \n(knots)', c='w')


    # Plot logo
    if (user_logo != 'n/a') & (user_logo[-4:] == '.png'):
        ax_logo = plt.axes([0, -0.12, 0.2, 0.2])
        logo = np.asarray(Image.open(user_logo))
        ax_logo.imshow(logo)
        ax_logo.set_xticks([])  # Remove x-axis ticks
        ax_logo.set_yticks([])  # Remove y-axis ticks
        ax_logo.set_facecolor('k')
        # ax_logo.set_aspect(1)
        # for spine in ax_logo.spines.values():
        #         spine.set_edgecolor('k')
    else:
        ax_logo = plt.axes([0, -0.12, 0.2, 0.2])
        ax_logo.set_xticks([])  # Remove x-axis ticks
        ax_logo.set_yticks([])  # Remove y-axis ticks
        ax_logo.set_facecolor('k')
        ax_logo.set_aspect(1)
        for spine in ax_logo.spines.values():
                spine.set_edgecolor('k')
        c1 = 'w'
        c2 = 'w'
        c3 = 'darkorange'
        lw1 = 3
        lw2 = 6
        ax_logo.plot([10,10], [10,25], linewidth=lw1, c=c1)
        ax_logo.plot([10,17], [10,10], linewidth=lw1, c=c1)
        ax_logo.plot([17,17], [10,25], linewidth=lw1, c=c1)

        ax_logo.plot([27,27], [10,25], linewidth=lw1, c=c1)
        ax_logo.plot([24,30], [25,25], linewidth=lw1, c=c1)
        ax_logo.plot([24,30], [10,10], linewidth=lw1, c=c1)

        ax_logo.plot([37,40], [17.5,17.5], linewidth=lw1, c=c1)

        ax_logo.plot([47,54], [25,25], linewidth=lw1, c=c1)
        ax_logo.plot([47,47], [25,17.5], linewidth=lw1, c=c1)
        ax_logo.plot([47,54], [17.5,17.5], linewidth=lw1, c=c1)
        ax_logo.plot([54,54], [10,17.5], linewidth=lw1, c=c1)
        ax_logo.plot([47,54], [10,10], linewidth=lw1, c=c1)

        ax_logo.plot([61,61], [10,25], linewidth=lw1, c=c1)
        ax_logo.plot([68,68], [10,25], linewidth=lw1, c=c1)
        ax_logo.plot([61,68], [10,10], linewidth=lw1, c=c1)
        ax_logo.plot([61,68], [25,25], linewidth=lw1, c=c1)

        ax_logo.plot([75,75], [10,25], linewidth=lw1, c=c1)
        ax_logo.plot([75,82], [25,10], linewidth=lw1, c=c1)
        ax_logo.plot([82,82], [10,25], linewidth=lw1, c=c1)

        ax_logo.plot([89,89], [10,25], linewidth=lw1, c=c1)
        ax_logo.plot([89,96], [25,22], linewidth=lw1, c=c1)
        ax_logo.plot([89,96], [10,13], linewidth=lw1, c=c1)
        ax_logo.plot([96,96], [13,22], linewidth=lw1, c=c1)

        ax_logo.plot([103,103], [10,25], linewidth=lw1, c=c1)
        ax_logo.plot([103,103], [10,25], linewidth=lw1, c=c1)
        ax_logo.plot([103,110], [10,10], linewidth=lw1, c=c1)
        ax_logo.plot([103,108], [17.5,17.5], linewidth=lw1, c=c1)
        ax_logo.plot([103,110], [25,25], linewidth=lw1, c=c1)

        # I
        ax_logo.plot([48,48], [52,93], linewidth=lw2, c=c2)
        ax_logo.plot([72,72], [52,93], linewidth=lw2, c=c2)
        ax_logo.plot([38,48], [52,52], linewidth=lw2, c=c2)
        ax_logo.plot([72,82], [52,52], linewidth=lw2, c=c2)
        ax_logo.plot([38,48], [93,93], linewidth=lw2, c=c2)
        ax_logo.plot([72,82], [93,93], linewidth=lw2, c=c2)
        ax_logo.plot([38,82], [108,108], linewidth=lw2, c=c2)
        ax_logo.plot([38,82], [37,37], linewidth=lw2, c=c2)
        ax_logo.plot([38,38], [93,108], linewidth=lw2, c=c2)
        ax_logo.plot([38,38], [37,52], linewidth=lw2, c=c2)
        ax_logo.plot([82,82], [93,108], linewidth=lw2, c=c2)
        ax_logo.plot([82,82], [37,52], linewidth=lw2, c=c2)

        ax_logo.fill_between([37, 82], [108, 108], [93, 93], color=c3)
        ax_logo.fill_between([37, 82], [52, 52], [38, 38], color=c3)
        ax_logo.fill_between([48, 72], [93, 93], [52, 52], color=c3)

        ax_logo.set_xlim([0,120])
        ax_logo.set_ylim([0,120])
        


    ## Plot thermodynamic parameters ##
    ax_params1 = plt.axes([0.22, -0.12, 0.58, 0.2])
    ax_params1.set_xticks([])  # Remove x-axis ticks
    ax_params1.set_yticks([])  # Remove y-axis ticks
    ax_params1.set_facecolor('k')
    for spine in ax_params1.spines.values():
            spine.set_edgecolor('w')

    # Titles
    ax_params1.text(0.07+0.04, 0.93, 'PCL', c='w', fontweight='bold', fontsize=14,
                    ha='center', va='top')
    ax_params1.text(0.23+0.04, 0.93, 'CAPE', c='w', fontweight='bold', fontsize=14,
                    ha='center', va='top')
    ax_params1.text(0.39+0.04, 0.93, 'CIN', c='w', fontweight='bold', fontsize=14,
                    ha='center', va='top')
    ax_params1.text(0.55+0.04, 0.93, 'LCL', c='w', fontweight='bold', fontsize=14,
                    ha='center', va='top')
    ax_params1.text(0.71+0.04, 0.93, 'LFC', c='w', fontweight='bold', fontsize=14,
                    ha='center', va='top')
    ax_params1.text(0.87+0.04, 0.93, 'EL', c='w', fontweight='bold', fontsize=14,
                    ha='center', va='top')
    ax_params1.set_xlim([0,1])
    ax_params1.set_ylim([0,1])
    ax_params1.plot([0, 1], [0.81, 0.81], c='w')
    # SFC
    if prof_choice == 'Surface':
        level_colors = ['red', 'lightskyblue', 'limegreen', 'yellow', 'magenta']
    else:
        level_colors = ['w', 'w', 'w', 'w', 'w']
    ax_params1.text(0.07+0.04, 0.75, 'SFC', c='w', fontsize=13.5,
                    ha='center', va='top')
    ax_params1.text(0.23+0.04, 0.75, f'{int(sfc_cape.magnitude)}', c=level_colors[0], fontsize=13.5,
                    ha='center', va='top')
    ax_params1.text(0.39+0.04, 0.75, f'{int(sfc_cin.magnitude)}', c=level_colors[1], fontsize=13.5,
                    ha='center', va='top')
    ax_params1.text(0.55+0.04, 0.75, f'{int(lcl_h.magnitude)}', c=level_colors[2], fontsize=13.5,
                    ha='center', va='top')
    if math.isnan(lfc_h.magnitude):
        ax_params1.text(0.71+0.04, 0.75, f'---', c=level_colors[3], fontsize=13.5,
                    ha='center', va='top')
    else:
        ax_params1.text(0.71+0.04, 0.75, f'{int(lfc_h.magnitude)}', c=level_colors[3], fontsize=13.5,
                        ha='center', va='top')
    if math.isnan(el_h.magnitude):
        ax_params1.text(0.87+0.04, 0.75, f'---', c=level_colors[4], fontsize=13.5,
                    ha='center', va='top')
    else:
        ax_params1.text(0.87+0.04, 0.75, f'{int(el_h.magnitude)}', c=level_colors[4], fontsize=13.5,
                        ha='center', va='top')
    # ML
    if prof_choice == 'Mixed Layer':
        level_colors = ['red', 'lightskyblue', 'limegreen', 'yellow', 'magenta']
    else:
        level_colors = ['w', 'w', 'w', 'w', 'w']
    ax_params1.text(0.07+0.04, 0.58, 'ML', c='w', fontsize=13.5,
                    ha='center', va='top')
    ax_params1.text(0.23+0.04, 0.58, f'{int(ml_cape.magnitude)}', c=level_colors[0], fontsize=13.5,
                    ha='center', va='top')
    ax_params1.text(0.39+0.04, 0.58, f'{int(ml_cin.magnitude)}', c=level_colors[1], fontsize=13.5,
                    ha='center', va='top')
    ax_params1.text(0.55+0.04, 0.58, f'{int(ml_lcl_h.magnitude)}', c=level_colors[2], fontsize=13.5,
                    ha='center', va='top')
    if math.isnan(ml_lfc_h.magnitude):
        ax_params1.text(0.71+0.04, 0.58, f'---', c=level_colors[3], fontsize=13.5,
                    ha='center', va='top')
    else:
        ax_params1.text(0.71+0.04, 0.58, f'{int(ml_lfc_h.magnitude)}', c=level_colors[3], fontsize=13.5,
                        ha='center', va='top')
    if math.isnan(ml_el_h.magnitude):
        ax_params1.text(0.87+0.04, 0.58, f'---', c=level_colors[4], fontsize=13.5,
                    ha='center', va='top')
    else:
        ax_params1.text(0.87+0.04, 0.58, f'{int(ml_el_h.magnitude)}', c=level_colors[4], fontsize=13.5,
                        ha='center', va='top')
    # MU
    if prof_choice == 'Most Unstable':
        level_colors = ['red', 'lightskyblue', 'limegreen', 'yellow', 'magenta']
    else:
        level_colors = ['w', 'w', 'w', 'w', 'w']
    ax_params1.text(0.07+0.04, 0.41, 'MU', c='w', fontsize=13.5,
                    ha='center', va='top')
    ax_params1.text(0.23+0.04, 0.41, f'{int(mu_cape.magnitude)}', c=level_colors[0], fontsize=13.5,
                    ha='center', va='top')
    ax_params1.text(0.39+0.04, 0.41, f'{int(mu_cin.magnitude)}', c=level_colors[1], fontsize=13.5,
                    ha='center', va='top')
    ax_params1.text(0.55+0.04, 0.41, f'{int(mu_lcl_h.magnitude)}', c=level_colors[2], fontsize=13.5,
                    ha='center', va='top')
    if math.isnan(mu_lfc_h.magnitude):
        ax_params1.text(0.71+0.04, 0.41, f'---', c=level_colors[3], fontsize=13.5,
                    ha='center', va='top')
    else:
        ax_params1.text(0.71+0.04, 0.41, f'{int(mu_lfc_h.magnitude)}', c=level_colors[3], fontsize=13.5,
                        ha='center', va='top')
    if math.isnan(mu_el_h.magnitude):
        ax_params1.text(0.87+0.04, 0.41, f'---', c=level_colors[4], fontsize=13.5,
                    ha='center', va='top')
    else:
        ax_params1.text(0.87+0.04, 0.41, f'{int(mu_el_h.magnitude)}', c=level_colors[4], fontsize=13.5,
                        ha='center', va='top')
    ax_params1.plot([0, 1], [0.25, 0.25], c='w')
    # Other thermo params
    ax_params1.text(0.05, 0.18, f'PWAT = {np.around(pw.magnitude, decimals=2)} in.', c='w', fontsize=13.5,
                    ha='left', va='top')
    if math.isnan(dcape.magnitude):
        ax_params1.text(0.37, 0.18, f'DCAPE = --- J/kg', c='w', fontsize=13.5,
                        ha='left', va='top')
    else:
        ax_params1.text(0.37, 0.18, f'DCAPE = {int(dcape.magnitude)} J/kg', c='w', fontsize=13.5,
                        ha='left', va='top')
#     if math.isnan(ship):
#         ax_params1.text(0.37, 0.18, f'SHIP = ---', c='w', fontsize=13.5,
#                         ha='left', va='top')
#     else:
#         ax_params1.text(0.37, 0.18, f'SHIP = {np.around(ship, decimals=1)}', c='w', fontsize=13.5,
#                         ha='left', va='top')
    

    ## Plot kinematic parameters ##
    ax_params2 = plt.axes([0.82, -0.12, 0.57, 0.2])
    ax_params2.set_xticks([])  # Remove x-axis ticks
    ax_params2.set_yticks([])  # Remove y-axis ticks
    ax_params2.set_facecolor('k')
    for spine in ax_params2.spines.values():
            spine.set_edgecolor('w')
    # Titles        
    ax_params2.text(0.32, 0.835, f'SRH (m$^2$/s$^2$)', c='w', fontweight='bold', fontsize=12,
                    ha='center', va='bottom')
    ax_params2.text(0.57, 0.85, f'Shear (kts)', c='w', fontweight='bold', fontsize=12,
                    ha='center', va='bottom')
    ax_params2.set_xlim([0,1])
    ax_params2.set_ylim([0,1])
    ax_params2.plot([0, 1], [0.81, 0.81], c='w')
    # 0-1
    ax_params2.text(0.015, 0.74, 'SFC-1km', c='w', fontsize=12.2,
                    ha='left', va='top')
    if math.isnan(srh01.magnitude):
        ax_params2.text(0.32, 0.74, f'---', c='w', fontsize=12.2,
                        ha='center', va='top')
    else:
        ax_params2.text(0.32, 0.74, f'{int(srh01.magnitude)}', c='w', fontsize=12.2,
                        ha='center', va='top')
    if math.isnan(shear01.magnitude):
        ax_params2.text(0.57, 0.74, f'---', c='w', fontsize=12.2,
                        ha='center', va='top')
    else:
        ax_params2.text(0.57, 0.74, f'{int(shear01.magnitude)}', c='w', fontsize=12.2,
                        ha='center', va='top')
    # 0-3
    ax_params2.text(0.015, 0.61, 'SFC-3km', c='w', fontsize=12.2,
                    ha='left', va='top')
    if math.isnan(srh03.magnitude):
        ax_params2.text(0.32, 0.61, f'---', c='w', fontsize=12.2,
                        ha='center', va='top')
    else:
        ax_params2.text(0.32, 0.61, f'{int(srh03.magnitude)}', c='w', fontsize=12.2,
                        ha='center', va='top')
    if math.isnan(shear03.magnitude):
        ax_params2.text(0.57, 0.61, f'---', c='w', fontsize=12.2,
                        ha='center', va='top')
    else:
        ax_params2.text(0.57, 0.61, f'{int(shear03.magnitude)}', c='w', fontsize=12.2,
                        ha='center', va='top')
    # EFF
    ax_params2.text(0.015, 0.48, 'Effect Inflow', c='w', fontsize=12.2,
                    ha='left', va='top')
    if math.isnan(srheff.magnitude):
        ax_params2.text(0.32, 0.48, f'---', c='w', fontsize=12.2,
                        ha='center', va='top')
    else:
        ax_params2.text(0.32, 0.48, f'{int(srheff.magnitude)}', c='cyan', fontsize=12.2,
                        ha='center', va='top')
    if math.isnan(sheareff.magnitude):
        ax_params2.text(0.57, 0.48, f'---', c='w', fontsize=12.2,
                        ha='center', va='top')
    else:
        ax_params2.text(0.57, 0.48, f'{int(sheareff.magnitude)}', c='cyan', fontsize=12.2,
                        ha='center', va='top')
    # 0-6
    ax_params2.text(0.015, 0.35, 'SFC-6km', c='w', fontsize=12.2,
                    ha='left', va='top')
    if math.isnan(shear06.magnitude):
        ax_params2.text(0.57, 0.35, f'---', c='w', fontsize=12.2,
                        ha='center', va='top')
    else:
        ax_params2.text(0.57, 0.35, f'{int(shear06.magnitude)}', c='w', fontsize=12.2,
                        ha='center', va='top')
    # 5-10
    ax_params2.text(0.015, 0.22, '5km-10km', c='w', fontsize=12.2,
                    ha='left', va='top')
    if math.isnan(shear510.magnitude):
        ax_params2.text(0.57, 0.22, f'---', c='w', fontsize=12.2,
                    ha='center', va='top')
    else:
        ax_params2.text(0.57, 0.22, f'{int(shear510.magnitude)}', c='w', fontsize=12.2,
                    ha='center', va='top')

    ax_params2.plot([0.73, 0.73], [0.81, 0], c='w')
    ax_params2.text(0.75, 0.7, 'Supercell:', c='w', fontsize=12.2, fontweight='bold',
                    ha='left', va='top')
    if math.isnan(supercell.magnitude):
        ax_params2.text(0.865, 0.55, f'---', c='w', fontsize=14, 
                    fontweight='bold', ha='center', va='top')
    else:
        if np.around(supercell.magnitude, decimals=1)[0] > 10:
            ax_params2.text(0.865, 0.55, f'{np.around(supercell.magnitude, decimals=1)[0]}', c='r', fontsize=14, 
                            fontweight='bold', ha='center', va='top')
        elif np.around(supercell.magnitude, decimals=1)[0] > 1:
            ax_params2.text(0.865, 0.55, f'{np.around(supercell.magnitude, decimals=1)[0]}', c='yellow', fontsize=14, 
                            fontweight='bold', ha='center', va='top')
        else:
            ax_params2.text(0.865, 0.55, f'{np.around(supercell.magnitude, decimals=1)[0]}', c='w', fontsize=14, 
                            fontweight='bold', ha='center', va='top')
    ax_params2.text(0.75, 0.35, 'STP (fixed):', c='w', fontsize=12.2, fontweight='bold',
                    ha='left', va='top')
    if math.isnan(stp.magnitude):
        ax_params2.text(0.865, 0.2, f'---', c='w', fontsize=14, 
                    fontweight='bold', ha='center', va='top')
    else:
        if np.around(stp.magnitude, decimals=1)[0] > 10:
            ax_params2.text(0.865, 0.2, f'{np.around(stp.magnitude, decimals=1)[0]}', c='magenta', fontsize=14, 
                            fontweight='bold', ha='center', va='top')
        elif np.around(stp.magnitude, decimals=1)[0] > 1:
            ax_params2.text(0.865, 0.2, f'{np.around(stp.magnitude, decimals=1)[0]}', c='r', fontsize=14, 
                            fontweight='bold', ha='center', va='top')
        elif np.around(stp.magnitude, decimals=1)[0] > 0.5:
            ax_params2.text(0.865, 0.2, f'{np.around(stp.magnitude, decimals=1)[0]}', c='yellow', fontsize=14, 
                            fontweight='bold', ha='center', va='top')
        else:
            ax_params2.text(0.865, 0.2, f'{np.around(stp.magnitude, decimals=1)[0]}', c='w', fontsize=14, 
                            fontweight='bold', ha='center', va='top')


    ## Plot hazard type
    ax_haz = plt.axes([1.41, -0.12, 0.18, 0.2])
    ax_haz.set_xticks([])  # Remove x-axis ticks
    ax_haz.set_yticks([])  # Remove y-axis ticks
    ax_haz.set_facecolor('k')
    ax_haz.set_xlim([0,1])
    ax_haz.set_ylim([0,1])
    for spine in ax_haz.spines.values():
            spine.set_edgecolor('w')
    ax_haz.plot([0, 1], [0.75, 0.75], c='w')
    ax_haz.text(0.5, 0.88, 'Psbl Haz. Type', c='w', fontweight='bold', fontsize=13.2,
                    ha='center', va='center')

    # Determine hazard
    if ((not math.isnan(supercell.magnitude)) & (not math.isnan(stp.magnitude))):
        if ((stp >= 3*units.dimensionless) & (srh01 >= 200*units('m^2/s^2')) & (srheff >= 200*units('m^2/s^2')) &
           (shear08 > 45*units.knots) & (lcl_h < 1000*units.m) & (ml_lcl_h < 1200*units.meters) &
           (LR01 >= 5) & (ml_cin > -50*units('J/kg')) & (hbot <= hgt[0])):
            hazard = 'PDS\nTOR'
            haz_color = 'magenta'
            # Missing SRW and Effective STP
        elif ((stp >= 4*units.dimensionless) & (ml_cin > -50*units('J/kg')) & (hbot <= hgt[0])):
            hazard = 'TOR'
            haz_color = 'r'
        elif ((stp >= 1*units.dimensionless) & (shear08 >= 40*units.knots) & (ml_cin > -50*units('J/kg')) &
             (hbot <= hgt[0])):
            hazard = 'TOR'
            haz_color = 'r'
            # Missing SRW option
        elif ((stp >= 1*units.dimensionless) & (LR01 >= 5) & (ml_cin > -50*units('J/kg')) & 
             (hbot <= hgt[0])):
            hazard = 'TOR'
            haz_color = 'r'
            # Missing Low and Mid RH avg
        elif ((stp >= 1*units.dimensionless) & (ml_cin > -150*units('J/kg')) & (hbot <= hgt[0])):
            hazard = 'MRGL\nTOR'
            haz_color = 'r'
        elif ((stp >= 0.5*units.dimensionless) & (ml_cin > -150*units('J/kg')) & (srh01 >= 150*units('m^2/s^2')) &
              (hbot <= hgt[0])):
            hazard = 'MRGL\nTOR'
            haz_color = 'r'
        elif (((stp >= 1*units.dimensionless) | (supercell >= 4*units.dimensionless)) & (mu_cin > -50*units('J/kg'))):
            hazard = 'SVR'
            haz_color = 'yellow'
        elif ((((supercell >= 2*units.dimensionless) & (ship >= 1)) | (dcape >= 750 * units('J/kg')) ) & 
              (mu_cin > -50*units('J/kg'))):
            hazard = 'SVR'
            haz_color = 'yellow'
#         elif ():
#             hazard = 'SVR'
#             haz_color = 'yellow'
            # Missing SigSvr and MMP
        elif ( ((wndg >= 0.5) | (ship >= 0.5) | (supercell >= 0.5*units.dimensionless)) & (mu_cin >= -75*units('J/kg')) ):
            hazard = 'MRGL\nSVR'
            haz_color = 'lightskyblue'
        elif (srf_HI > 105 * units.degF):
            hazard = 'EXCESSIVE\nHEAT'
            haz_color = 'orange'
        # MORE...
        else:
            hazard = 'NONE'
            haz_color = 'palegoldenrod'
    else:
        hazard = 'NONE'
        haz_color = 'palegoldenrod'

    ax_haz.text(0.5, 0.4, hazard, c=haz_color, fontweight='bold', fontsize=17,
                    ha='center', va='center')

    #fig.tight_layout()
    fig.patch.set_facecolor('k')

    # Save file to same location as original data file with identical name
    #savefile = 'Desktop/UI_Sonde/SkewTs/'+year+month+day+'_'+hour+minute+'_SkewT.png'
    savefile = filepath.split('.txt')[0] + '.png'
    plt.savefig(f'{savefile}', dpi=dpi_choice, bbox_inches='tight')
    plt.close()

    os.system(f"open {savefile}")



def import_file():
    global sounding_file_path
    file_path = filedialog.askopenfilename(title="Select a file", filetypes=[("Text files", "*.txt"), ("All files", "*.*")])
    if file_path:
        #print("Selected file:", file_path)
        sounding_file_path = file_path
        
logo_file_path = 'n/a'
def import_logo():
    global logo_file_path
    file_path2 = filedialog.askopenfilename(title="Select a file", filetypes=[("PNGs", "*.png"), ("All files", "*.*")])
    if file_path2:
        #print("Selected file:", file_path)
        logo_file_path = file_path2


def submit_callback(e1, e2, sounding_file_path, combobox, dpi, logo_file_path):
#     print("User entered : " + e1.get())
#     print("User entered : " + e2.get())
    #print("Will now run file with : " + e1.get() + e2.get() + sounding_file_path)
#     os.system(f'python Desktop/UI_Sonde/UI_Sonde_run.py --inputlat={e1.get()} --inputlon={e2.get()} --inputfile={sounding_file_path}')
    plot_sonde(e1.get(), e2.get(), sounding_file_path, combobox.get(), dpi, logo_file_path)
    return None


# Creating tkinter window
window = tk.Tk()
window.title('UI-Sonde')
#window.geometry('450x670')

input_label = tk.Label(window, text="Welcome to UI-Sonde!",
                       font = ("Ariel", 20, 'bold'), foreground="navy").grid(column = 0,row = 1, padx = 0, pady = 25, 
                                                                     columnspan = 2)
description_label = tk.Label(window, text="Select a file and input launch site lat/lon to generate a Skew-T",
                             foreground="black", font = ("Ariel", 12)).grid(column = 0,row = 2, padx = 10, 
                                                                            pady = 25, columnspan = 2)

# Add Logo
# img = PhotoImage(file = 'Desktop/UI_Sonde/CliMAS_Logo.png')
# img1 = img.subsample(2, 2)
# Label(window, image = img1).grid(row = 1, column = 1,
#        columnspan = 2, padx = 5, pady = 5)


# Input values
Label(window, text='Station Latitude (degrees)').grid(row=5, padx=15, sticky='E')
Label(window, text='Station Longitude (degrees)').grid(row=6, padx=15, sticky='E')
e1 = Entry(window)
e2 = Entry(window)
e1.grid(row=5, column=1, padx=10)
e2.grid(row=6, column=1, padx=10)


submit_button = tk.Button(window, text = "Generate!" , font = ("Ariel", 13, 'bold', 'underline'), foreground = 'blue',
                          command = lambda: submit_callback(e1, e2,   sounding_file_path, combobox, dpi, 
                                                            logo_file_path))
submit_button.grid(column=1, row=7, pady=5)

# Create an "Import File" button
import_button = tk.Button(window, text="Import Sonde", command=import_file)
import_button.grid(row=4, column=1, pady=20)

description_label = tk.Label(window, text="Optional Settings...\n___________________________________________________________",
                             foreground="black", font = ("Ariel", 12, 'italic')).grid(column = 0,row = 8, padx = 20, pady = 15, columnspan=2)

# Drop down for profile type
Label(window, text='Parcel Type').grid(row=9, padx=15, pady=15, sticky='E')
combobox= ttk.Combobox(window, state= "readonly")
combobox['values']=('Surface','Mixed Layer','Most Unstable')
combobox.current(0)
combobox.grid(row=9, column=1, pady=15)

# Select dpi
dpi = 100
def get_choice(returned_val):
    global dpi
    dpi = returned_val

Label(window, text='Image Resolution').grid(row=10, padx=15, sticky='E')
Label(window, text='Default: Medium (100 dpi)', font = ("Ariel", 8)).grid(row=11, padx=15, sticky='E')
#v = StringVar(window, "Medium") 
v = IntVar()
v.set(100)
r1 = Radiobutton(window, text="Low", value=50, var=v, 
                 command=lambda *args: get_choice(50)).grid(row=10, column=1, sticky='W')
r2 = Radiobutton(window, text="Medium", value=100, var=v, 
                 command=lambda *args: get_choice(100)).grid(row=11, column=1, sticky='W')
r3 = Radiobutton(window, text="High", value=200, var=v, 
                 command=lambda *args: get_choice(200)).grid(row=12, column=1, sticky='W')

# # Choose radar on/off
# radar = 0
# def get_choice2(returned_val2):
#     global radar
#     radar = returned_val2

# Label(window, text='').grid(row=13, padx=15, pady=2, sticky='E')
# Label(window, text='Plot Nearest Radar').grid(row=14, padx=15, sticky='E')
# Label(window, text='Internet Connection Required', font = ("Ariel", 8)).grid(row=15, padx=15, sticky='E')
# v2 = IntVar()
# v2.set(0)
# c1 = Radiobutton(window, text="On", var=v2, value=1,
#                  command=lambda *args: get_choice2(1)).grid(row=14, column=1, sticky='W')
# c2 = Radiobutton(window, text="Off", var=v2, value=0,
#                  command=lambda *args: get_choice2(0)).grid(row=15, column=1, sticky='W')


# Create an "Import Logo" button
Label(window, text='Import Custom Logo').grid(row=16, padx=15, sticky='E')
import_button = tk.Button(window, text="Select File", command=import_logo)
import_button.grid(row=16, column=1, pady=15, sticky='W')


tk.Label(window, text="UI-Sonde v0.4.0",
                             foreground="orange", font = ("Ariel", 12)).grid(column = 0,row = 17, padx = 10, 
                                                                            pady = 20, columnspan = 2)

def open_new_window():
    new_window = Toplevel(window)  # Create a new window
    new_window.title("About")
    #new_window.geometry("300x300")  

    Label(new_window, text="About UI-Sonde",
                       font = ("Ariel", 20, 'bold'), foreground="navy", justify='center').grid(column = 0,row = 1, 
                                                                                               padx = 20, pady = 25, 
                                                                                               columnspan = 5)
    text1 = ("UI-Sonde is a lightweight MacOS and Windows app for plotting Skew-Ts in the SHARPpy format."+
            " UI-Sonde accepts the same radiosonde file format as SHARPpy and is capable of generating"+
            " PNG Skew-T plots with a similar GUI interface. The software is built entirely in Python,"+
            " utlizing the MetPy package for most calculations. Additional calculations were adapted"+
            " from SHARPpy or built specifically for UI-Sonde.")
    Label(new_window, text=text1, font = ("Ariel", 12), justify='center', wraplength=250).grid(pady=20, column=0, row=2,
                                                                                           columnspan=5, padx=20)
    
#about = tk.Label(window, text="About", cursor='hand1', foreground='navy').grid(row=15, pady=10)
#about.bind("<Button-1>", lambda e: open_new_window())
tk.Button(window, text="About", font = ("Ariel", 10, 'italic'), command=open_new_window, width=5, cursor='hand1',
          fg='black').grid(row=18, column=0, columnspan=2)


mainloop()