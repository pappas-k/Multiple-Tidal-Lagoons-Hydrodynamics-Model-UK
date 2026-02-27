import uptide
import uptide.tidal_netcdf
import datetime
from math import tanh
import sys                                  #af
sys.path.append('../')                      #af
from inputs.simulation_parameters import *  #af
from tools import utm
utm_zone = i_zone
utm_band = i_band

# def initial_forcing(t_start):
constituents = ['Q1', 'O1', 'P1', 'K1', 'N2', 'M2', 'S2', 'K2']
tide = uptide.Tides(constituents)
tide.set_initial_time(datetime.datetime(s_year,s_month,s_day,s_hour,s_min))
tnci = uptide.tidal_netcdf.OTPSncTidalInterpolator(tide,grid_forcing_file,hf_forcing_file, 
                                                    ranges=range_forcing_coords)

def set_tidal_field(elev, t, start):
  tnci.set_time(t)
  mesh2d = elev.function_space().mesh()
  xvector = mesh2d.coordinates.dat.data
  evector = elev.dat.data
  ramp = tanh((t-start)/25000.)
  for i,xy in enumerate(xvector):
    lat, lon = utm.to_latlon(xy[0], xy[1], utm_zone, utm_band)
    try:
      evector[i] = tnci.get_val((lon, lat))*ramp    # Adding initial a correction depth for LAT
    except uptide.netcdf_reader.CoordinateError:
      evector[i] = 0.    # Adding initial a correction depth for LAT
