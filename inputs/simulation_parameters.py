"""
Input file paths
"""
# TODO: This is a good start, but perhaps a json file format or similar might be more generic going forward
#model_data_dir = '/home/chien/Documents/model_data'
model_data_dir ="../../model_data"
#model_data_dir = '/Home/Dropbox/0_PhD_projects/model_data'

# Mesh
#mesh_file = 'inputs/ambient_UK_mesh10000.msh'
mesh_file = "inputs/west_uk_multiple_lagoons.msh"  # Elevation boundaries - 1,2
#mesh_file = "inputs/west_uk_multiple_lagoons_coarser.msh"  # Elevation boundaries - 1,2



# Bathymetry file(s) - list in order of highest -> lowest resolution, second entry is source of data
# CHECK PREPROCESSING SCRIPT ON ASSUMPTIONS ABOUT THESE BATHYMETRIES
#bathymetry_file= f"{model_data_dir}/cropped_irish-sea.nc" # source is z [digimap]

#bathymetry_file_0 = f"{model_data_dir}/digimap_West_UK_6_arc_seconds.nc" # source is z [digimap]
bathymetry_file_0 = f"{model_data_dir}/digimap_West_UK_1_arc_seconds.nc" # source is z [digimap]
bathymetry_file_1 = f"{model_data_dir}/GEBCO_West_UK.nc" # source is elevation [emod]

# bathymetry_file_0 = f"{model_data_dir}/cropped_irish-sea.nc" # source is z [digimap]
# bathymetry_file_1 = f"{model_data_dir}/emod_data/emod_west_UK.nc" # source is z [emod]

# TODO: Suggested format to consider, on allocating info for bathymetry
# bathymetries = [(bathymetry_file, source_variable (e.g. 'z','Band1'...), Datum check (e.g. 'MWL'...)]
bathymetries = [(bathymetry_file_0, 'Band1', 'LAT'),
                (bathymetry_file_1, 'elevation', 'MWL')]

# Forcing
grid_forcing_file = f'{model_data_dir}/gridES2008.nc'
hf_forcing_file = f'{model_data_dir}/hf.ES2008.nc'
range_forcing_coords = ((-12., -2.), (48, 59))

# Detectors
i_tidegauge_file = 'inputs/useful_gauges_BODC.csv'
elevation_detectors = []
additional_detector_files = ['inputs/extra_detectors_TRS']

# Bed morphology data file
bed_classification_file = f'{model_data_dir}/BGS_data/bed_class_pentland_rev.nc'

# Friction file (e.g. for variable friction)
friction_data = "inputs/n_max_125.npy"
use_friction_data = False

"""
Outputs folders
"""
# Outputs folder
ramp_output_folder = 'outputs/outputs_ramp'
run_output_folder = 'outputs/outputs_run'

"""
UTM parameters - update based on location considered
"""
i_zone = 30   # used in detectors.py, tidal_amplitude_serial.py, bathymetry.py
i_band = 'U'  # used in tidal_amplitude_serial.py, bathymetry.py - 30N is northern hemisphere, V describes North Sea

"""
Simulation start time parameters (for TPXO)
"""
s_year = 2002
s_month = 10
s_day = 20
s_hour = 0
s_min = 0

"""
Detector parameters
"""
# Maximum distance from detectors
max_dist = 5e3

"""
Eikonal and preprocessing stage parameters 
"""
# Characteristic length for eikonal eqt
i_L = 1e3

# epss values set the accuracy (in meters) of the final "distance to boundary" function.
i_epss = [100000., 10000., 5000., 2500., 1500., 1000.]

# Boundary values (from QGIS)
open_bnd = [4,5,6]
internal_lagoon_bnd = [16,17,26,27,36,37,46,47,56,57,66,67,76,77]
external_lagoon_bnd = [18,19,28,29,38,39,48,49,58,59,68,69,78,79]
land_bnd = 1000
inner_id = [12,22,32,42,52,62,72]
outer_id = [13,23,33,43,53,63,73]

#names of lagoons
lagoons = ['SW', 'CA', 'WA', 'CO', 'LI', 'BL', 'SO']

#hydraulic structures boundaries ids
lagoon_bnd={
            'SW': {'internal_tb_bnd' : 16, 'internal_sl_bnd': 17, 'external_tb_bnd' : 18 , 'external_sl_bnd' : 19, 'inner_id' : 12, 'outer_id' : 13 },
            'CA': {'internal_tb_bnd' : 26, 'internal_sl_bnd': 27, 'external_tb_bnd' : 28 , 'external_sl_bnd' : 29, 'inner_id' : 22, 'outer_id' : 23 },
            'WA': {'internal_tb_bnd' : 36, 'internal_sl_bnd': 37, 'external_tb_bnd' : 38 , 'external_sl_bnd' : 39, 'inner_id' : 32, 'outer_id' : 33 },
            'CO': {'internal_tb_bnd' : 46, 'internal_sl_bnd': 47, 'external_tb_bnd' : 48 , 'external_sl_bnd' : 49, 'inner_id' : 42, 'outer_id' : 43 },
            'LI': {'internal_tb_bnd' : 56, 'internal_sl_bnd': 57, 'external_tb_bnd' : 58 , 'external_sl_bnd' : 59, 'inner_id' : 52, 'outer_id' : 53 },
            'BL': {'internal_tb_bnd' : 66, 'internal_sl_bnd': 67, 'external_tb_bnd' : 68 , 'external_sl_bnd' : 69, 'inner_id' : 62, 'outer_id' : 63 },
            'SO': {'internal_tb_bnd' : 76, 'internal_sl_bnd': 77, 'external_tb_bnd' : 78 , 'external_sl_bnd' : 79, 'inner_id' : 72, 'outer_id' : 73 },
            }

"""
Lagoon parameters
"""
turbine_params={'f_g': 50, 'g_p': 96, 'g': 9.807, 't_d': 7.35, 't_cap': 20,
                'h_cap': 5, 'dens': 1025, 'h_min': 1.0, 'eta': [0.95, 0.95],
                'options': 1}
N_t = 30.
N_s = 15.
operation_ramp = 'two-way'
operation = 'two-way'
#operation = "two-way-pump"



"""
Bathymetry parameters
"""
# Minimum depth in bathymetry
i_min_depth = -10.0     # used in bathymetry.py, tidal_amplitude_serial.py
# Name of bathymetry parameter in the nc file
alt_name = 'z'

"""
Ramp and run parameters
"""
i_ramptime = 2 * 24 * 3600      # Ramptime in seconds
i_dt = 100                      # Crank-Nicolson timestep in seconds
i_alpha = 1.5                   # Wetting and drying parameter
i_manning = 0.03                # Manning parameter
i_lat_cor = 53                  # Coriolis calculation parameters, latitude in degrees

i_spatial_harmonics_distribution = False
ramp_exp_interval = 1000.      # Ramp output interval
run_exp_interval = 10000.         # Run output interval
i_t_end = 30 * 24 * 3600         # Simulation duration in seconds
# i_t_end = 3 * 24 * 3600         # Test on Docker

"""
Other parameters (based on specific thetis application)
"""
# Shear stresses parameters (in run.py)
grav_acc = 9.807
density  = 1025

# tide constituents, sorted by amplitude
i_constituents = ['M2', 'S2', 'N2', 'K2', 'K1', 'O1', 'P1', 'Q1']
# i_constituents = ['Q1', 'O1', 'P1', 'K1', 'N2', 'M2', 'S2', 'K2']  # unsorted
