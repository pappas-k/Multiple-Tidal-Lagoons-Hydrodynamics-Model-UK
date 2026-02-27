from tools import utm
from netCDF4 import Dataset as NetCDFFile
import scipy.interpolate
import numpy as np
from firedrake import *

minimum_depth = -10.0

utm_zone=30
utm_band='N'


def get_val(X,t):
  lat, lon = utm.to_latlon(X[0], X[1], utm_zone, utm_band)
  return max(minimum_depth, -nci.get_val((lat, lon)))


# TODO: Merged iterative and original get_bathymetry script - We might need to think about handling particular datasets
def get_bathymetry(bathymetry_file, mesh2d, source='z', bathymetry_function=None):
    """
    Function to produce bathymetry
    :param bathymetry_file: Bathymetry file
    :param mesh2d: domain discretisation as a Firedrake mesh object
    :param source: Variable name used for interpolation quantity
    :param bathymetry_function: Use of pre-existing bathymetry function, in case this is an iterative process
    :return:
    """
    nc = NetCDFFile(bathymetry_file)
    print(bathymetry_file)
    print(nc.variables.keys())
    lat = np.float64(nc.variables['lat'][:])
    lon = np.float64(nc.variables['lon'][:])
    values = np.float64(nc.variables[source][:, :])
    values = values.filled(9999.)
    interpolator = scipy.interpolate.RegularGridInterpolator((lat, lon), values)
    P1_2d = FunctionSpace(mesh2d, 'CG', 1)
    if bathymetry_function is None:
        bathymetry_function = Function(P1_2d, name="bathymetry").assign(np.nan)
    xvector = mesh2d.coordinates.dat.data
    bvector = bathymetry_function.dat.data
    assert xvector.shape[0] == bvector.shape[0]
    for a, xy in enumerate(xvector):
        lat, lon = utm.to_latlon(xy[0], xy[1], utm_zone, utm_band)
        try:
            bvector[a] = max(-interpolator((lat, lon)), minimum_depth)
            # if -interpolator((lat, lon)) > 2 * minimum_depth and np.isnan(bvector[a]):
            #     bvector[a] = -interpolator((lat, lon))
        except ValueError:
            bvector[a] = np.nan
    return bathymetry_function


def get_bed_class(class_file, mesh2d, default=6):
    nc = NetCDFFile(class_file)
    lat = nc.variables['lat'][:]
    lon = nc.variables['lon'][:]
    values = nc.variables['Band1'][:,:]
    values = values.filled(9999.)
    interpolator = scipy.interpolate.RegularGridInterpolator((lat, lon), values)
    P1_2d = FunctionSpace(mesh2d, 'CG', 1)
    bedclass2d = Function(P1_2d, name="bed_class")
    xvector = mesh2d.coordinates.dat.data
    bvector = bedclass2d.dat.data
    assert xvector.shape[0]==bvector.shape[0]
    for i,xy in enumerate(xvector):
        lat, lon = utm.to_latlon(xy[0], xy[1], utm_zone, utm_band)
        try:
            bvector[i] = interpolator((lat, lon))
        except ValueError:
            bvector[i] = default
    return bedclass2d


def get_manning_class(class_file, mesh2d, manning_function, default=6,):
    nc = NetCDFFile(class_file)
    lat = nc.variables['lat'][:]
    lon = nc.variables['lon'][:]
    values = nc.variables['Band1'][:,:]
    values = values.filled(9999.)
    interpolator = scipy.interpolate.RegularGridInterpolator((lat, lon), values,  method='nearest')
    P1_2d = FunctionSpace(mesh2d, 'CG', 1)
    manningclass2d = Function(P1_2d, name="manning")
    xvector = mesh2d.coordinates.dat.data
    bvector = manningclass2d.dat.data
    assert xvector.shape[0]==bvector.shape[0]
    for i,xy in enumerate(xvector):
        lat, lon = utm.to_latlon(xy[0], xy[1], utm_zone, utm_band)
        try:
            if interpolator((lat,lon)) >= 1:
                bvector[i] = manning_function(interpolator((lat, lon))-1)
            elif interpolator((lat,lon)) < 1:
                bvector[i] = manning_function(default)
        except ValueError:
            bvector[i] = manning_function(default)
    return manningclass2d


def get_bathymetry_iteration(bathymetry_file, mesh2d, bathymetry_function, source='digimap'):
    print("Iterating bathymetry")
    nc = NetCDFFile(bathymetry_file)
    lat = nc.variables['lat'][:]
    lon = nc.variables['lon'][:]
    if source=='digimap':
        values = nc.variables['z'][:,:]
    else:
        values = nc.variables['elevation'][:,:]
    values = values.filled(9999.)
    interpolator = scipy.interpolate.RegularGridInterpolator((lat, lon), values)
    xvector = mesh2d.coordinates.dat.data
    bvector = bathymetry_function.dat.data
    print(np.isnan(np.sum(bvector)))
    assert xvector.shape[0]==bvector.shape[0]
    for i,xy in enumerate(xvector):
        lat, lon = utm.to_latlon(xy[0], xy[1], utm_zone, utm_band)

        if np.isnan(bvector[i]):
            print("replacing")
            try:
                bvector[i] = max(-interpolator((lat, lon)), minimum_depth)
            except ValueError:
                bvector[i] = np.nan
    return bathymetry_function


def smoothen_bathymetry(bathymetry2d):
  v = TestFunction(bathymetry2d.function_space())
  massb = assemble(v * bathymetry2d *dx)
  massl = assemble(v*dx)
  with massl.dat.vec as ml, massb.dat.vec as mb, bathymetry2d.dat.vec as sb:
      ml.reciprocal()
      sb.pointwiseMult(ml, mb)


def get_bathymetry_from_text(bath, txtfile):
    mesh2d = bath.function_space().mesh()
    xvector = mesh2d.coordinates.dat.data
    evector = bath.dat.data
    data = np.loadtxt(txtfile)
    intp = scipy.interpolate.NearestNDInterpolator(data[:, 0:2], data[:, 2])
    for i, xy in enumerate(xvector):
        evector[i] = intp(xy)
    return bath
