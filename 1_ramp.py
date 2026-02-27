# 2D shallow water equations
# ================================================================================
# A ramp simulation of the multiple lagoons model
# ================================================================================
from thetis import *
import tools.tidal_forcing_ramp
import tools.thetis_support_scripts
import sys  # af

sys.path.append('../')  # af
from inputs.simulation_parameters import *  # af
from firedrake.petsc import PETSc  # af
from modules import input_barrages
from modules.tools import LagoonCallback
import numpy as np
import warnings  # af

warnings.simplefilter(action="ignore", category=DeprecationWarning)  # af

inputdir = 'inputs'
outputdir = ramp_output_folder
datadir = 'data'

with timed_stage('reading mesh'):
    mesh2d = Mesh(mesh_file)

PETSc.Sys.Print('Loaded mesh ' + mesh2d.name)  # af
PETSc.Sys.Print('Exporting to ' + outputdir)  # af

# simulation ID
identifier = -1
PETSc.Sys.Print('Simulation identifier : ' + str(identifier))  # af

ramptime = i_ramptime
t_start = - ramptime  # Simulation start time relative to tidal_forcing
t_end = ramptime + t_start  # Simulation duration in sec
Dt = i_dt  # Time integrator timestep
t_export = ramp_exp_interval  # Export time if necessary
wd_alpha = i_alpha  # Wetting and drying
mu_manning = Constant(i_manning)  # Bed friction

lat_coriolis = i_lat_cor  # Coriolis calculation parameters
CG_2d = FunctionSpace(mesh2d, 'CG', 1)
bathymetry_2d, h_viscosity = tools.thetis_support_scripts.initialise_fields(mesh2d, inputdir,
                                                                            outputdir, identifier,
                                                                            manning=False)
# mu_manning = Constant(0.026)
coriolis_2d = tools.thetis_support_scripts.coriolis(mesh2d, lat_coriolis)

# Lagoon initialisation # lagoon_control, lagoon_params
lagoon_input = input_barrages.input_predefined_barrage_specs(N_t, N_s, operation=operation_ramp,
                                                             turbine_parameters=turbine_params)

with timed_stage('initialisation'):
    # --- create solver ---
    solver_obj = solver2d.FlowSolver2d(mesh2d, bathymetry_2d)
    options = solver_obj.options
    options.cfl_2d = 1.0
    options.use_nonlinear_equations = True
    options.simulation_export_time = t_export
    options.simulation_end_time = ramptime
    options.coriolis_frequency = coriolis_2d
    options.output_directory = outputdir
    options.check_volume_conservation_2d = True
    options.fields_to_export = ['elev_2d', 'uv_2d']
    options.fields_to_export_hdf5 = []
    options.element_family = "dg-dg"
    options.swe_timestepper_type = 'CrankNicolson'
    options.swe_timestepper_options.implicitness_theta = 0.6
    options.swe_timestepper_options.use_semi_implicit_linearization = True
    options.use_wetting_and_drying = True
    options.wetting_and_drying_alpha = Constant(wd_alpha)
    options.manning_drag_coefficient = mu_manning
    options.horizontal_viscosity = h_viscosity
    options.use_grad_div_viscosity_term = True
    options.use_grad_depth_viscosity_term = False
    options.timestep = Dt  # override dt for CrankNicolson (semi-implicit)
    options.swe_timestepper_options.solver_parameters = {
        'snes_type': 'newtonls',
        'snes_rtol': 1e-2,  # 1e-3
        'snes_linesearch_type': 'bt',
        'snes_max_it': 20,  # 20
        'ksp_type': 'preonly',
        'pc_type': 'lu',
        'pc_factor_mat_solver_package': 'mumps',
    }

# Initialising hydraulic structures
lagoon_hydraulic_structures = {'SW': {"sl_i": Constant(0.), "sl_o": Constant(0.), "tb_i": Constant(0.), "tb_o": Constant(0.)},
                                'CA': {"sl_i": Constant(0.), "sl_o": Constant(0.), "tb_i": Constant(0.), "tb_o": Constant(0.)},
                                'WA': {"sl_i": Constant(0.), "sl_o": Constant(0.), "tb_i": Constant(0.), "tb_o": Constant(0.)},
                                'CO': {"sl_i": Constant(0.), "sl_o": Constant(0.), "tb_i": Constant(0.), "tb_o": Constant(0.)},
                                'LI': {"sl_i": Constant(0.), "sl_o": Constant(0.), "tb_i": Constant(0.), "tb_o": Constant(0.)},
                                'BL': {"sl_i": Constant(0.), "sl_o": Constant(0.), "tb_i": Constant(0.), "tb_o": Constant(0.)},
                                'SO': {"sl_i": Constant(0.), "sl_o": Constant(0.), "tb_i": Constant(0.), "tb_o": Constant(0.)}
                                 }

tidal_elev = Function(bathymetry_2d.function_space())
solver_obj.bnd_functions['shallow_water'] = {
                                            open_bnd[0]: {'elev': tidal_elev},
                                            open_bnd[1]: {'elev': tidal_elev},
                                            open_bnd[2]: {'elev': tidal_elev}
                                            }

for lag in lagoons:
    solver_obj.bnd_functions['shallow_water'][lagoon_bnd[lag]['internal_tb_bnd']] = {'flux': lagoon_hydraulic_structures[lag]["tb_i"]}
    solver_obj.bnd_functions['shallow_water'][lagoon_bnd[lag]['external_tb_bnd']] = {'flux': lagoon_hydraulic_structures[lag]["tb_o"]}
    solver_obj.bnd_functions['shallow_water'][lagoon_bnd[lag]['internal_sl_bnd']] = {'flux': lagoon_hydraulic_structures[lag]["sl_i"]}
    solver_obj.bnd_functions['shallow_water'][lagoon_bnd[lag]['external_sl_bnd']] = {'flux': lagoon_hydraulic_structures[lag]["sl_o"]}


# Create a dat file for lagoon operation outputs
f1 = open(datadir + "/Lagoon_" + str(identifier) + ".dat", "w")

elev_init = Function(CG_2d).assign(0.0)

solver_obj.assign_initial_conditions(uv=as_vector((1e-3, 0.0)), elev=elev_init)

for lag in lagoons:
    cb_lagoon = LagoonCallback(solver_obj, {"inner": dx(lagoon_bnd[lag]['inner_id']), "outer": dx(lagoon_bnd[lag]['outer_id'])}, lagoon_input,
                               thetis_boundaries=lagoon_hydraulic_structures[lag], time=t_start,
                               number_timesteps=5, name="lagoon_" + str(lag) + "_" + str(identifier),   # number_timesteps=5
                               adaptive_control=None, export_to_hdf5=True)  #delete the outputs_ramp files to avoid error:  File "h5py/h5f.pyx", line 126, in h5py.h5f.create
                                                                                                                           #BlockingIOError: [Errno 11] Unable to create file (unable to lock file, errno = 11, error message = 'Resource temporarily unavailable')
                                                                                                                           #application called MPI_Abort(PYOP2_COMM_WORLD, 1) - process 0

    solver_obj.add_callback(cb_lagoon, 'timestep')

uv, elev = solver_obj.timestepper.solution.split()

# P1_2d = FunctionSpace(mesh2d, 'CG', 1)

def intermediate_steps(t):
    # Temporary output to data file.
    # for item in cb_lagoon.output:
    #     f1.write(" {:8.3f} ".format(item))
    # f1.write("\n")
    # f1.flush()

    # Exporting to data file - useful for quick sampling etc.
    if i_spatial_harmonics_distribution == True \
            and t % ramp_exp_interval == 0:
        PETSc.Sys.Print("Exporting elevation field for harmonic analysis")  # af
        elev_CG = Function(CG_2d, name='elev_CG').project(elev)
        checkpoint_file = checkpointing.DumbCheckpoint(outputdir + '/elev_' + str(t))
        checkpoint_file.store(elev_CG)
        checkpoint_file.close()

    # Export final state that can be picked up later - like load state but including lagoon state
    if t == t_end: tools.thetis_support_scripts.export_final_state(inputdir, identifier, uv, elev,
                                                                   lagoon=[cb_lagoon.status])


def update_forcings(t):
    intermediate_steps(float(t + t_start))
    PETSc.Sys.Print("Updating tidal field at t={}".format(t_start + t))  # af
    tools.tidal_forcing_ramp.set_tidal_field(tidal_elev, t + int(t_start), t_start)


solver_obj.iterate(update_forcings=update_forcings)
