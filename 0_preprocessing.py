"""
Pre-processing script. This script:
- Interpolates bathymetry to an hdf5 file that can be imported later
- Adds viscosity sponges in boundary conditions
- Can be used to edit bathymetry/ friction and other fields that then feed into the simulations

In modifying the regions close to boundary the Eikonal equation is solved
using Firedrake (initial script from Roan, modified by Stephan and Than)
"""

from scipy.interpolate import interp1d
from thetis import *
from tools import bathymetry, tidal_amplitude
import inputs.simulation_parameters as inputs
import numpy as np

import sys
from firedrake.petsc import PETSc           #af



inputdir = "inputs"
outputdir = "outputs"
mesh = Mesh(inputs.mesh_file)

# Step 0 - Calculate lowest astronomical tide if bathymetry dataset is based on Lowest Astronomical Tide
# This applies for digimap datasets or digitised maps from admiralty sources
V = FunctionSpace(mesh, 'CG', 1)

lat = Function(V)
tidal_amplitude.get_lowest_astronomical_tide(lat)
File(outputdir + '/lat.pvd').write(lat)

# Step 1 - Calculate distance for viscosity
PETSc.Sys.Print("Calculate distance for viscosity") #af

# Boundary conditions
bcs = []
for i in inputs.open_bnd:
    bcs.append(DirichletBC(V, 0.0, i))
# for j in inputs.internal_lagoon_bnd:
#     bcs.append(DirichletBC(V, 1e5, j))

for lag in inputs.lagoons:
    x = inputs.lagoon_bnd[lag]['internal_tb_bnd']
    y = inputs.lagoon_bnd[lag]['internal_sl_bnd']
    bcs.append(DirichletBC(V, 1e5, x))
    bcs.append(DirichletBC(V, 1e5, y))



L = inputs.i_L
v = TestFunction(V)
u = Function(V)

solver_parameters = {
    'snes_type': 'ksponly',
    'ksp_rtol': 1e-4,
    'ksp_type': 'preonly',
    'pc_type': 'lu',
    'pc_factor_mat_solver_packages': 'mumps',
    }

# Before we solve the Eikonal equation, let's solve a Laplace equation to
# generate an initial guess
F = L**2*(inner(grad(u), grad(v))) * dx - v * dx
solve(F == 0, u, bcs, solver_parameters=solver_parameters)
solver_parameters = {
    'snes_type': 'newtonls',
    'ksp_rtol': 1e-4,
        'ksp_type': 'preonly',
        'pc_type': 'lu',
        'pc_factor_mat_solver_packages': 'mumps',
        }

epss = inputs.i_epss
for i, eps in enumerate(epss):
  PETSc.Sys.Print("Solving Eikonal with eps == ", float(eps))   #af
  F = inner(sqrt(inner(grad(u), grad(u))), v) * dx - v * dx + eps*inner(grad(u), grad(v)) * dx
  solve(F == 0, u, bcs, solver_parameters=solver_parameters)

File(outputdir + "/dist.pvd").write(u)

"""
Adding viscosity sponge
"""
chk = DumbCheckpoint(inputdir + "/viscosity", mode=FILE_CREATE)
with timed_stage('initialising viscosity'):
    h_viscosity = Function(V, name="viscosity")
    h_viscosity.interpolate(max_value(1., 1000 * (1. - u / 2e4)))
    chk.store(h_viscosity, name="viscosity")
    File(outputdir + '/viscosity.pvd').write(h_viscosity)

"""
Creating a Manning/Quadratic drag/other type friction field to be used in the simulations
"""
with timed_stage('initialising manning friction'):
    if inputs.use_friction_data:
        manning_data = np.load(inputs.friction_data)
        interpolator = np.vectorize(interp1d(manning_data[0, :], manning_data[1, :],
                                             fill_value=(manning_data[1, 0], manning_data[1, -1]),
                                             bounds_error=False))
        manning_2d = bathymetry.get_manning_class(inputs.bed_classification_file, mesh, interpolator)
        bathymetry.smoothen_bathymetry(manning_2d)
        File(outputdir + '/manning.pvd').write(manning_2d)
        print_output('Exported manning')
    else:
        manning_2d = Function(V, name='manning')
        manning_2d.interpolate(max_value(inputs.i_manning, 0.1 * (1. - u / 5e4)))
        File(outputdir + '/manning.pvd').write(manning_2d)

        # # Apply Manning in the specified rectangular patch
        # x, y = SpatialCoordinate(mesh)
        # #coordinates of the patch
        # xmin = 284294
        # xmax = 328812
        # ymin = 5771112
        # ymax = 5792432
        # new_manning=  0.1
        # #new_manning = Constant(0.1)
        # #manning_patch = Function(V).interpolate(new_manning)
        # manning_export = Function(V)
        # manning_export.interpolate(conditional(And(And(ge(x, xmin), le(x, xmax)), And(ge(y, ymin), le(y, ymax))),
        #                                        new_manning, manning_2d))
        #
        # File(outputdir + '/manning.pvd').write(manning_export)

"""
Interpolating bathymetry
"""
# Note in LAT datum cases, the 'lat' field is added to correct to mean water level
# This assumes that Gebco is always the least fine bathymetry and is corrected for MWL already and that all Meygen &
# Digimap files need to be corrected
with timed_stage('initialising bathymetry'):
    bath=None
    for i, (f, source, datum) in enumerate(inputs.bathymetries):
        bath = bathymetry.get_bathymetry(f, mesh, source=source, bathymetry_function=bath)
        if datum=='LAT':
            bath.assign(bath + lat)
        File(outputdir + '/bath'+str(i)+'.pvd').write(bath)


    def zbedf(bathy, distance):
        """
        Function used to edit bathymetry close to a boundary (e.g. when wetting and drying needs to be avoided)
        :param bathy:  Bathymetry field
        :param distance:  Distance from particular boundary determined by the Eikonal equation
        :return: Modified bathymetry field
        """
        zbed = conditional(ge(bathy,25 *(1.- distance / 15000.)),bathy,
                           25 *(1.- distance / 15000.))
        return zbed

    # Applying bathymetry correction at the boundary
    #
    #bath.interpolate(max_value(zbedf(bath,u), -100.))
    bath.interpolate(max_value(zbedf(bath, u), inputs.i_min_depth))


    L = 1e3
    # Set up hydraulic structures
    bc2 = []
    #for i in inputs.internal_lagoon_bnd:
    #    bc2.append(DirichletBC(V, 0.0, i))
    #for j in inputs.external_lagoon_bnd:
    #    bc2.append(DirichletBC(V, 0.0, j))

    for lag in inputs.lagoons:
        i = inputs.lagoon_bnd[lag]['internal_tb_bnd']
        j = inputs.lagoon_bnd[lag]['internal_sl_bnd']
        k = inputs.lagoon_bnd[lag]['external_tb_bnd']
        l = inputs.lagoon_bnd[lag]['external_sl_bnd']
        bc2.append(DirichletBC(V, 0, i))
        bc2.append(DirichletBC(V, 0, j))
        bc2.append(DirichletBC(V, 0, k))
        bc2.append(DirichletBC(V, 0, l))

    v2 = TestFunction(V)
    u2 = Function(V)

    solver_parameters = {
        'snes_type': 'ksponly',
        'ksp_rtol': 1e-4,
        'ksp_type': 'preonly',
        'pc_type': 'lu',
        'pc_factor_mat_solver_packages': 'mumps',
    }
    # # Before we solve the Eikonal equation, let's solve a Laplace equation to
    # # generate an initial guess
    F = L ** 2 * (inner(grad(u2), grad(v2))) * dx - v2 * dx
    solve(F == 0, u2, bc2, solver_parameters=solver_parameters)
    solver_parameters = {
        'snes_type': 'newtonls',
        'ksp_rtol': 1e-4,
        'ksp_type': 'preonly',
        'pc_type': 'lu',
        'pc_factor_mat_solver_packages': 'mumps',
    }

    epss = inputs.i_epss
    for i, eps in enumerate(epss):
        print("Solving Eikonal with eps == ", float(eps))
        F = inner(sqrt(inner(grad(u2), grad(u2))), v2) * dx - v2 * dx + eps * inner(grad(u2), grad(v2)) * dx
        solve(F == 0, u2, bc2, solver_parameters=solver_parameters)


    # File("outputs/dist.pvd").write(u2)

    # Testing conditional statements   / bath has been introduced previously
    def zbedf(bathy, distance):
        zbed = conditional(le(distance, 100), 30.0 - distance / 100 * 5,
                           conditional(le(distance, 200), 25,
                                       conditional(le(distance, 1000), max_value(25 - (distance - 200) / 1000 * 40, bathy),
                                                   bathy)))
        return zbed


    bath.interpolate(zbedf(bath, u2))  # testing conditional statements
    # Smoothing bathymetry
    bathymetry.smoothen_bathymetry(bath)

    V_DG = FunctionSpace(mesh, 'DG', 1)
    bathDG=Function(V_DG).project(bath)
    File(outputdir + '/bath_DG.pvd').write(bathDG)
    File(outputdir + '/bath.pvd').write(bath)
    chk = DumbCheckpoint(inputdir + "/bathymetry2D", mode=FILE_CREATE)
    chk.store(bath, name="bathymetry")