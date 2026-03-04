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

import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.mpl.geoaxes as geo

from PIL import Image


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Inputs"
    )
    parser.add_argument("--inputlat", required=True, type=float)
    parser.add_argument("--inputlon", required=True, type=float)
    parser.add_argument("--inputfile", required=True, type=str)
    args = parser.parse_args()

    station_lat = args.inputlat
    station_lon = args.inputlon
    user_filepath = args.inputfile


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


# Parse input data

launch_lat = station_lat
launch_lon = station_lon

filepath = user_filepath


# Create map image to plot
fig = plt.figure(figsize=(8, 12))              
ax_map = plt.axes(projection=ccrs.PlateCarree())
ax_map.add_feature(cfeature.STATES, linewidth=3, edgecolor='w')
ax_map.add_feature(cfeature.LAND, facecolor='k')
ax_map.add_feature(cfeature.OCEAN)
#ax_map.add_feature(USCOUNTIES, linewidth=0.5, edgecolor='gray')
ax_map.set_extent([launch_lon-8, launch_lon+8, launch_lat-5, launch_lat+5])
ax_map.scatter(launch_lon, launch_lat, color='r', s=1500, marker='*')
plt.savefig('Temp_Map_Inset.png', bbox_inches='tight', dpi=100)
plt.close()

    

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
lcl_h = hgt_fit(p, hgt, lcl_pressure)
el_h = hgt_fit(p, hgt, el_pressure)
lfc_h =hgt_fit(p, hgt, lfc_pressure)
sfc_cape, sfc_cin = mpcalc.surface_based_cape_cin(p, T, Td)
# ML
ml_p, ml_T, ml_Td = mpcalc.mixed_parcel(p, T, Td)
ml_prof = mpcalc.parcel_profile(p, ml_T, ml_Td).to('degC')
ml_lcl_pressure, ml_lcl_temperature = mpcalc.lcl(ml_p, ml_T, ml_Td)
ml_el_pressure, _ = mpcalc.el(p, T, Td, parcel_temperature_profile = ml_prof)
ml_lfc_pressure, _ = mpcalc.lfc(p, T, Td, ml_prof, which='most_cape')
ml_lcl_h = hgt_fit(p, hgt, ml_lcl_pressure)
ml_el_h = hgt_fit(p, hgt, ml_el_pressure)
ml_lfc_h =hgt_fit(p, hgt, ml_lfc_pressure)
ml_cape, ml_cin = mpcalc.mixed_layer_cape_cin(p, T, Td)
# MU
mu_p, mu_T, mu_Td, mu_ind = mpcalc.most_unstable_parcel(p, T, Td)
mu_prof = mpcalc.parcel_profile(p[mu_ind:], mu_T, mu_Td).to('degC')
mu_lcl_pressure, mu_lcl_temperature = mpcalc.lcl(mu_p, mu_T, mu_Td)
mu_el_pressure, _ = mpcalc.el(p[mu_ind:], T[mu_ind:], Td[mu_ind:], parcel_temperature_profile = mu_prof)
mu_lfc_pressure, _ = mpcalc.lfc(p[mu_ind:], T[mu_ind:], Td[mu_ind:], mu_prof, which='most_cape')
mu_lcl_h = hgt_fit(p, hgt, mu_lcl_pressure)
mu_el_h = hgt_fit(p, hgt, mu_el_pressure)
mu_lfc_h =hgt_fit(p, hgt, mu_lfc_pressure)
mu_cape, mu_cin = mpcalc.most_unstable_cape_cin(p, T, Td)

# Calculate additional parameters
pw = mpcalc.precipitable_water(p, Td).to(units.inches)
#dcape, _, _ = mpcalc.downdraft_cape(p, T, Td)

if hgt[-1:][0] < 6000*units.meters:
    rm = [0 * units.knots, 0 * units.knots]
    lm = [0 * units.knots, 0 * units.knots]
    mean_storm = [0 * units.knots, 0 * units.knots]
else:
    rm, lm, mean_storm = mpcalc.bunkers_storm_motion(p, u, v, hgt)
# Shear and SRH
if hgt[-1:][0] < 1000*units.meters:
    shear01 = np.nan * units.knots
    srh01 = np.nan * units('m/s^2')
else:
    u01, v01 = mpcalc.bulk_shear(p, u, v, height = hgt, depth = 1000*units.meters, bottom = hgt[0])
    shear01 = np.sqrt((u01**2) + (v01**2))
    srh01,_,_ = mpcalc.storm_relative_helicity(hgt, u, v, depth=1000*units.meters, storm_u=rm[0], storm_v=rm[1])

if hgt[-1:][0] < 3000*units.meters:
    shear03 = np.nan * units.knots
    srh03 = np.nan * units('m/s^2')
else:
    u03, v03 = mpcalc.bulk_shear(p, u, v, height = hgt, depth = 3000*units.meters, bottom = hgt[0])
    shear03 = np.sqrt((u03**2) + (v03**2))
    srh03,_,_ = mpcalc.storm_relative_helicity(hgt, u, v, depth=3000*units.meters, storm_u=rm[0], storm_v=rm[1])

if hgt[-1:][0] < 6000*units.meters:
    shear06 = np.nan * units.knots
else:
    u06, v06 = mpcalc.bulk_shear(p, u, v, height = hgt, depth = 6000*units.meters, bottom = hgt[0])
    shear06 = np.sqrt((u06**2) + (v06**2))

if hgt[-1:][0] < 10000*units.meters:
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
if hgt[-1:][0] < 8000*units.meters:
    shear08 = np.nan * units.knots
else:
    u08, v08 = mpcalc.bulk_shear(p, u, v, height = hgt, depth = 8000*units.meters, bottom = hgt[0])
    shear08 = np.sqrt((u08**2) + (v08**2))
#
if hgt[-1:][0] < 3500*units.meters:
    wndg = np.nan * units.dimensionless
else:
    ind_3km = np.where(hgt>=3000*units.meters)[0][0]
    LR03 = -1/3 * ((T[ind_3km]-T[0]).magnitude)
    ind_1km = np.where(hgt>=1000*units.meters)[0][0]
    ind_3_5km = np.where(hgt>=3500*units.meters)[0][0]
    
    Mean_Wind = np.mean(mpcalc.wind_speed(u[ind_1km:ind_3_5km], v[ind_1km:ind_3_5km]).to(units('m/s')))
    wndg = (wind_damage(ml_cape, ml_cin, LR03, Mean_Wind).magnitude)*units.dimensionless
    

    
    

#### Plotting time! ####

# Make figure
fig = plt.figure(figsize=(9, 9), facecolor='k', constrained_layout=True)
#fig, ax = plt.subplots(figsize=(9, 9), constrained_layout=True) 

# Plot Skew-T
skew = SkewT(fig, rotation=45)
skew.plot(p, T, 'r', linewidth=2)
skew.plot(p, Td, 'limegreen', linewidth=2)
skew.plot_barbs(p[::25], u[::25], v[::25], barbcolor='w')
skew.ax.set_ylim(1000, 90)
skew.ax.set_xlim(-50, 50)

skew.ax.set_xlabel(f'Temperature ({T.units:~P})', fontsize=12, c='w')
skew.ax.set_ylabel(f'Pressure ({p.units:~P})', fontsize=12, c='w')
skew.ax.tick_params(labelsize=12, labelcolor='w', color='w')
skew.ax.set_facecolor('k')

# Plot parcel profile
skew.plot(p, prof, 'w', linewidth=2, linestyle='--')
#skew.plot(p[mu_ind:], mu_prof, 'w', linewidth=2, linestyle='--')

# Add lines for lcl, lfc, el, effective
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

#Shade areas of CAPE and CIN
skew.shade_cin(p, T, prof, Td)
skew.shade_cape(p, T, prof)

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
ax_hod.set_title('Plotted with UI-Sonde v1.0 + MetPy', loc='right', c='w')


# Plot map
ax_map_img = plt.axes([1.2, 0.06, 0.35, 0.35])
img = np.asarray(Image.open('Temp_Map_Inset.png'))
ax_map_img.imshow(img)
ax_map_img.set_xticks([])  # Remove x-axis ticks
ax_map_img.set_yticks([])  # Remove y-axis ticks
ax_map_img.set_facecolor('k')
ax_map_img.set_aspect(1)
for spine in ax_map_img.spines.values():
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
ax_logo = plt.axes([0, -0.12, 0.2, 0.2])
logo = np.asarray(Image.open('CliMAS_Logo.png'))
ax_logo.imshow(logo)
ax_logo.set_xticks([])  # Remove x-axis ticks
ax_logo.set_yticks([])  # Remove y-axis ticks
ax_logo.set_facecolor('k')
# ax_logo.set_aspect(1)
# for spine in ax_logo.spines.values():
#         spine.set_edgecolor('k')
        

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
ax_params1.text(0.07+0.04, 0.75, 'SFC', c='w', fontsize=13.5,
                ha='center', va='top')
ax_params1.text(0.23+0.04, 0.75, f'{int(sfc_cape.magnitude)}', c='w', fontsize=13.5,
                ha='center', va='top')
ax_params1.text(0.39+0.04, 0.75, f'{int(sfc_cin.magnitude)}', c='w', fontsize=13.5,
                ha='center', va='top')
ax_params1.text(0.55+0.04, 0.75, f'{int(lcl_h.magnitude)}', c='limegreen', fontsize=13.5,
                ha='center', va='top')
if math.isnan(lfc_h.magnitude):
    ax_params1.text(0.71+0.04, 0.75, f'---', c='yellow', fontsize=13.5,
                ha='center', va='top')
else:
    ax_params1.text(0.71+0.04, 0.75, f'{int(lfc_h.magnitude)}', c='yellow', fontsize=13.5,
                    ha='center', va='top')
if math.isnan(el_h.magnitude):
    ax_params1.text(0.87+0.04, 0.75, f'---', c='magenta', fontsize=13.5,
                ha='center', va='top')
else:
    ax_params1.text(0.87+0.04, 0.75, f'{int(el_h.magnitude)}', c='magenta', fontsize=13.5,
                    ha='center', va='top')
# ML
ax_params1.text(0.07+0.04, 0.58, 'ML', c='w', fontsize=13.5,
                ha='center', va='top')
ax_params1.text(0.23+0.04, 0.58, f'{int(ml_cape.magnitude)}', c='w', fontsize=13.5,
                ha='center', va='top')
ax_params1.text(0.39+0.04, 0.58, f'{int(ml_cin.magnitude)}', c='w', fontsize=13.5,
                ha='center', va='top')
ax_params1.text(0.55+0.04, 0.58, f'{int(ml_lcl_h.magnitude)}', c='w', fontsize=13.5,
                ha='center', va='top')
if math.isnan(lfc_h.magnitude):
    ax_params1.text(0.71+0.04, 0.58, f'---', c='w', fontsize=13.5,
                ha='center', va='top')
else:
    ax_params1.text(0.71+0.04, 0.58, f'{int(ml_lfc_h.magnitude)}', c='w', fontsize=13.5,
                    ha='center', va='top')
if math.isnan(el_h.magnitude):
    ax_params1.text(0.87+0.04, 0.58, f'---', c='w', fontsize=13.5,
                ha='center', va='top')
else:
    ax_params1.text(0.87+0.04, 0.58, f'{int(ml_el_h.magnitude)}', c='w', fontsize=13.5,
                    ha='center', va='top')
# MU
ax_params1.text(0.07+0.04, 0.41, 'MU', c='w', fontsize=13.5,
                ha='center', va='top')
ax_params1.text(0.23+0.04, 0.41, f'{int(mu_cape.magnitude)}', c='w', fontsize=13.5,
                ha='center', va='top')
ax_params1.text(0.39+0.04, 0.41, f'{int(mu_cin.magnitude)}', c='w', fontsize=13.5,
                ha='center', va='top')
ax_params1.text(0.55+0.04, 0.41, f'{int(mu_lcl_h.magnitude)}', c='w', fontsize=13.5,
                ha='center', va='top')
if math.isnan(lfc_h.magnitude):
    ax_params1.text(0.71+0.04, 0.41, f'---', c='w', fontsize=13.5,
                ha='center', va='top')
else:
    ax_params1.text(0.71+0.04, 0.41, f'{int(mu_lfc_h.magnitude)}', c='w', fontsize=13.5,
                    ha='center', va='top')
if math.isnan(el_h.magnitude):
    ax_params1.text(0.87+0.04, 0.41, f'---', c='w', fontsize=13.5,
                ha='center', va='top')
else:
    ax_params1.text(0.87+0.04, 0.41, f'{int(mu_el_h.magnitude)}', c='w', fontsize=13.5,
                    ha='center', va='top')
ax_params1.plot([0, 1], [0.25, 0.25], c='w')
# Other thermo params
ax_params1.text(0.05, 0.18, f'PWAT = {np.around(pw.magnitude, decimals=2)} in.', c='w', fontsize=13.5,
                ha='left', va='top')
# ax_params1.text(0.37, 0.18, f'DCAPE = {int(dcape.magnitude)} J/kg', c='w', fontsize=13.5,
#                 ha='left', va='top')


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
    elif ((supercell >= 2*units.dimensionless) & (mu_cin > -50*units('J/kg'))):
        hazard = 'SVR'
        haz_color = 'yellow'
        # Missing DCAPE and SHIP (can add fairly easily)
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

savefile = year+month+day+'_'+hour+minute+'_SkewT.png'
plt.savefig(f'{savefile}', dpi=100, bbox_inches='tight')
plt.close()

os.system(f"open {savefile}")