import numpy as np
import inputs.input_file_paths
import datetime
import uptide
import math
from firedrake import *
sys.path.append('../')                      #af
from tools.processing_support_scripts import output_field_h5
from tools.signal_processing import produce_reconstructed_signal, determine_mean_tidal_range,theoretical_tidal_stream_power
from netCDF4 import Dataset as NetCDFFile
from scipy.interpolate import RegularGridInterpolator

from mpi4py import MPI
from matplotlib import pyplot as plt


comm = MPI.COMM_WORLD

# inputdir = 'inputs'
# outputdir = inputs.input_file_paths.output_folder

l = int(sys.argv[1])
outputdir, inputdir, friction_factor = inputs.input_file_paths.cases(l=l)
# outputdir ='outputs'

mesh2d = Mesh(inputs.input_file_paths.mesh_file)
xvector = mesh2d.coordinates.dat.data

# Initialising functions and functionspaces
P1_2D = FunctionSpace(mesh2d, "CG", 1)
P1v_2D = VectorFunctionSpace(mesh2d, 'CG', 1)
elev = Function(P1_2D, name='elev_CG')
uv = Function(P1v_2D, name="uv_CG")
bathymetry_2d = Function(P1_2D, name="bathymetry")

M2_amp = Function(P1_2D, name='M2_amp')
S2_amp = Function(P1_2D, name='S2_amp')
M2_phase = Function(P1_2D, name='M2_phase')
S2_phase = Function(P1_2D, name='S2_phase')
range_field = Function(P1_2D, name='tidal_range')
velocity_field = Function(P1_2D, name='mean_velocity_magnitude')
stream_field = Function(P1_2D, name='tidal_stream_energy')
depth_field = Function(P1_2D, name='mean_depth')
range_energy_field = Function(P1_2D, name='tidal_range_energy')
stream_energy_field = Function(P1_2D, name='tidal_stream_energy')
stream_power_field = Function(P1_2D, name='tidal_stream_power')
stream_power_field_computed = Function(P1_2D, name='tidal_stream_power_computed')
diameter_field = Function(P1_2D, name='feasible_diameter')

# Idealised tidal stream turbine characteristics
diameter = 25
turbine_area = lambda d: np.pi * d**2 / 4

# Set simulation start time
dt = inputs.input_file_paths.elevation_output_interval
datetime_start = datetime.datetime(*inputs.input_file_paths.simulation_time)

# Set time periods of harmonic analysis
t_start = 0
# t_end = t_start + 29.5 * 24 * 3600
t_end = t_start + 29.5 * 24 * 3600
t_n = int((t_end - t_start)/ dt + 1)

thetis_times = t_start + dt * np.arange(t_n) + dt
constituents = ['M2', 'S2', 'Q1', 'O1', 'P1', 'K1', 'N2',  'K2']

# Reading bathymetry
chk = DumbCheckpoint(inputdir+"/bathymetry2D", mode=FILE_READ)
chk.load(bathymetry_2d)
chk.close()
bath_data_set = bathymetry_2d.dat.data

# Initialising arrays for elevation and velocities
print(outputdir + '/elev_' +str(int(t_start+dt)) + '.0.h5')
checkpoint_file = checkpointing.DumbCheckpoint(outputdir + '/elev_' +str(int(t_start+dt)) + '.0', mode=FILE_READ)
checkpoint_file.load(elev)
checkpoint_file.close()
elev_data_set = np.empty((t_n, elev.dat.data.shape[0]))
u_data_set = np.empty((t_n, elev.dat.data.shape[0]))
v_data_set = np.empty((t_n, elev.dat.data.shape[0]))
print(elev.dat.data.shape[0])

for i in range(t_n):
    print('Reading h5 files. Time  ', i, len(range(t_n)))
    if comm.rank == 0: print('Reading Elevation')
    chk = checkpointing.DumbCheckpoint(outputdir+'/elev_' + str(t_start+(i + 1) * dt), mode=FILE_READ)
    chk.load(elev)
    chk.close()
    elev_data_set[i, :] = elev.dat.data[:]

    if comm.rank == 0: print('Reading Velocities')
    chk = DumbCheckpoint(outputdir + "/uv_" + str(t_start+(i + 1) * dt), mode=FILE_READ)
    chk.load(uv)
    chk.close()
    u_data_set[i,:] = uv.dat.data[:,0]
    v_data_set[i,:] = uv.dat.data[:,1]


detector_amplitudes = []
detector_phases = []
tidal_range, range_energy = [], []
mean_depth, mean_stream_power, mean_stream_power_computed, stream_energy, mean_velocity = [], [], [], [], []

for i in range(elev.dat.data.shape[0]):
    thetis_elev = elev_data_set[:, i]
    thetis_u = u_data_set[:, i]
    thetis_v = v_data_set[:, i]
    tide = uptide.Tides(constituents)
    tide.set_initial_time(datetime.datetime(*inputs.input_file_paths.simulation_time))

    predicted_time_elevation_series = produce_reconstructed_signal(thetis_times[:], thetis_elev[:], constituents)
    range, potential_energy = determine_mean_tidal_range(predicted_time_elevation_series[0],
                                                         predicted_time_elevation_series[1])
    predicted_time_u_series = produce_reconstructed_signal(thetis_times[:], thetis_u[:],constituents)
    predicted_time_v_series = produce_reconstructed_signal(thetis_times[:], thetis_v[:],constituents)
    predicted_signal_uv_mag = np.sqrt(predicted_time_u_series[1]**2+predicted_time_v_series[1]**2)
    computed_magnitude = np.sqrt(thetis_u[:]**2+thetis_v[:]**2)

    stream_power_computed = theoretical_tidal_stream_power(computed_magnitude[:], diameter=diameter)/(turbine_area(diameter)*1e3)
    stream_power = theoretical_tidal_stream_power(predicted_signal_uv_mag, diameter=diameter)/(turbine_area(diameter) * 1e3)
    average_stream_power = np.mean(stream_power)
    average_stream_power_computed =np.mean(stream_power_computed)

    if comm.rank == 0: print(i, average_stream_power_computed, average_stream_power)

    # Subtract mean
    thetis_elev = thetis_elev - thetis_elev.mean()
    thetis_amplitudes, thetis_phases = uptide.analysis.harmonic_analysis(tide, thetis_elev[:], thetis_times[:])

    mean_velocity.append(float(np.average(computed_magnitude)))
    mean_stream_power.append(float(average_stream_power))
    mean_stream_power_computed.append((float(average_stream_power_computed)))
    mean_depth.append(float(bath_data_set[i] - thetis_elev.mean()))
    detector_amplitudes.append(thetis_amplitudes)
    detector_phases.append(thetis_phases)
    tidal_range.append(float(range))
    range_energy.append(float(potential_energy/1e6))


range_energy_field.dat.data[:] = np.array(range_energy)[:]
range_field.dat.data[:] = np.array(tidal_range)[:]
depth_field.dat.data[:] = np.array(mean_depth)[:]
velocity_field.dat.data[:] = np.array(mean_velocity)[:]
stream_power_field.dat.data[:] = np.array(mean_stream_power)[:]
stream_power_field_computed.dat.data[:] = np.array(mean_stream_power_computed)[:]

M2_amp.dat.data[:] = np.array(detector_amplitudes)[:, 0]
M2_phase.dat.data[:] = np.array(detector_phases)[:, 0]

S2_amp.dat.data[:] = np.array(detector_amplitudes)[:, 1]
S2_phase.dat.data[:] = np.array(detector_phases)[:, 1]

File('outputs/range.pvd').write(range_field, range_energy_field, depth_field)
File('outputs/stream.pvd').write(velocity_field, stream_power_field, stream_power_field_computed)
File('outputs/amp.pvd').write(M2_amp, S2_amp)
File('outputs/phase.pvd').write(M2_phase, S2_phase)
# M2_phase.dat.data[:] = np.arcsin(np.sin(M2_phase.dat.data[:]))
M2_phase.dat.data[:] = np.remainder(M2_phase.dat.data[:],2*math.pi)*360/(2*math.pi)
# S2_phase.dat.data[:] = np.arcsin(np.sin(S2_phase.dat.data[:]))
S2_phase.dat.data[:] = np.remainder(S2_phase.dat.data[:],2*math.pi)*360/(2*math.pi)
File('outputs/phase_mod_pi.pvd').write(M2_phase, S2_phase)

output_field_h5(outputdir,range_field,"tidal_range")
output_field_h5(outputdir,range_energy_field, "tidal_range_energy")
output_field_h5(outputdir,M2_amp,"M2_amp")
output_field_h5(outputdir,M2_phase,"M2_phase")
output_field_h5(outputdir,S2_amp,"S2_amp")
output_field_h5(outputdir,S2_phase,"S2_phase")

#
# import utm
#
# # Read TPXO dataset to compare to
# hRe = NetCDFFile('netcdf/hf.ES2008.nc').variables['hRe'][:]
# hIm = NetCDFFile('netcdf/hf.ES2008.nc').variables['hIm'][:]
# lon = NetCDFFile('netcdf/hf.ES2008.nc').variables['lon_z'][:]
# lat = NetCDFFile('netcdf/hf.ES2008.nc').variables['lat_z'][:]
#
# TPXO_M2_amplitude = np.sqrt(np.square(hRe) + np.square(hIm))
# TPXO_M2_amp_interp = RegularGridInterpolator((lon[:, 0], lat[0, :]), TPXO_M2_amplitude[0, :, :])
#
# # Interpolate onto mesh
# TPXO_M2 = Function(P1_2D)
# TPXO_M2_data = TPXO_M2.dat.data
# for i, xy in enumerate(xvector):
#     lat, lon = utm.to_latlon(xy[0], xy[1], 30, 'U', strict=False)
#     try:
#         TPXO_M2_data[i] = TPXO_M2_amp_interp((lon, lat))
#     except ValueError:
#         TPXO_M2_data[i] = 0.0
#
# File('tpxo_amp.pvd').write(TPXO_M2)
#
# abs_diff = Function(P1_2D)
# abs_diff.dat.data[:] = TPXO_M2_data - M2_amp.dat.data[:]
#
# File('absolute_difference.pvd').write(abs_diff)
#
# rel_diff = Function(P1_2D)
# rel_diff.dat.data[:] = (TPXO_M2_data - M2_amp.dat.data[:]) / np.maximum(TPXO_M2_data, 0.01)
#
# File('relative_difference.pvd').write(rel_diff)