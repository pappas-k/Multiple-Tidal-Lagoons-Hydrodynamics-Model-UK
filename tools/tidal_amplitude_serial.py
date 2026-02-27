"""
This script takes the expected tidal amplitude across the domain
and outputs it for a transformation between datums.
"""
import uptide
import uptide.tidal_netcdf
import datetime
# import utm
import numpy
from netCDF4 import Dataset as NetCDFFile
from firedrake import *
import scipy.interpolate
import os
import sys                                  #af
sys.path.append('../')                      #af
from inputs.simulation_parameters import *  #af
from tools import utm

os.chdir('../')

utm_zone = i_zone
utm_band = i_band
minimum_depth = i_min_depth

#def initial_forcing(t_start):
constituents = ['Q1', 'O1', 'P1', 'K1', 'N2', 'M2', 'S2', 'K2']
tide = uptide.Tides(constituents)
tide.set_initial_time(datetime.datetime(s_year,s_month,s_day,s_hour,s_min))
tnci = uptide.tidal_netcdf.OTPSncTidalInterpolator(tide, grid_forcing_file, hf_forcing_file,
                                                    ranges=range_forcing_coords)

def write_lowest_astronomical_tide(mesh2d):
   amp = numpy.sqrt(tnci.real_part**2 + tnci.imag_part**2)
   val = amp.sum(axis=0)
   tnci.interpolator = uptide.netcdf_reader.Interpolator(tnci.nci.origin, tnci.nci.delta, val, tnci.nci.mask)
   xvector = mesh2d.coordinates.dat.data
   data = []
   for i,xy in enumerate(xvector):
     lat, lon = utm.to_latlon(xy[0], xy[1], utm_zone, utm_band)
     try:
       val = tnci.get_val((lon, lat), allow_extrapolation=True)
       data.append([xy[0], xy[1], val])
     except uptide.netcdf_reader.CoordinateError:
       pass
   numpy.savetxt('inputs/lat.txt', data)


def get_bathymetry_serial(bathymetry_file, mesh2d):
   nc = NetCDFFile(bathymetry_file)
   lat = nc.variables['lat'][:]
   lon = nc.variables['lon'][:]
   values = nc.variables[alt_name][:,:]
   values = values.filled(9999.)
   interpolator = scipy.interpolate.RegularGridInterpolator((lat, lon), values)
   P1_2d = FunctionSpace(mesh2d, 'CG', 1)
   bathymetry2d = Function(P1_2d, name="bathymetry")
   xvector = mesh2d.coordinates.dat.data
   bvector = bathymetry2d.dat.data
   data = []
   assert xvector.shape[0]==bvector.shape[0]
   for i,xy in enumerate(xvector):
       lat, lon = utm.to_latlon(xy[0], xy[1], utm_zone, utm_band)
       bvector[i] = max(-interpolator((lat, lon)), minimum_depth)
       data.append([xy[0], xy[1], bvector[i]])
   numpy.savetxt('inputs/bathy.txt', data)


"""
Make sure you pick the right file for the bathymetry below -
"""
mesh = Mesh(mesh_file)
write_lowest_astronomical_tide(mesh)
get_bathymetry_serial(bathymetry_file, mesh)
