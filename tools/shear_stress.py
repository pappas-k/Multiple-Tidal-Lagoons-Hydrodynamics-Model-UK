# TODO this is work in progress
from thetis import *
import sys                                                              #af
sys.path.append('../')                                                  #af
from inputs.simulation_parameters import *                              #af


def wd_bathymetry_displacement(solver,functionspace):
    """
    Returns wetting and drying bathymetry displacement as described in:
    Karna et al.,  2011.
    """
    H = solver.fields["bathymetry_2d"]+solver.fields["elev_2d"]
    disp = Function(functionspace).assign(0.5 * (sqrt(H ** 2 + solver.options.wetting_and_drying_alpha ** 2) - H))
    return disp


def compute_total_depth(solver,functionspace):
    """
    Returns effective depth by accounting for the wetting and drying algorithm
    """
    if hasattr(solver.options, 'use_wetting_and_drying') and solver.options.use_wetting_and_drying:
        return Function(functionspace).assign(solver.fields["bathymetry_2d"]+solver.fields["elev_2d"]+
                                              wd_bathymetry_displacement(solver, functionspace))
    else:
        return Function(functionspace).assign(solver.fields["bathymetry_2d"]+solver.fields["elev_2d"])


def compute_bed_shear_term(solver,functionspace):
    # Parameters
    g_grav = grav_acc
    dens = density
    uv, elev = solver.timestepper.solution.split()

    C_D = 0
    if solver_obj.options.quadratic_drag_coefficient is not None:
        C_D = Function(functionspace).assign(solver_obj.options.quadratic_drag_coefficient)
    elif solver_obj.options.manning_drag_coefficient is not None:
        C_D = Function(functionspace).assign(g_grav* solver_obj.options.manning_drag_coefficient**2/
                                             compute_total_depth(solver,functionspace)**(1./3.))

    shear_stress_vector = Function(VectorFunctionSpace(solver.mesh2d, "DG",1)).\
        interpolate(conditional(le(elev+solver.fields["bathymetry_2d"],0), as_vector((0.0,0.0)),dens * C_D * sqrt(dot(uv,uv)) * uv))
    return shear_stress_vector