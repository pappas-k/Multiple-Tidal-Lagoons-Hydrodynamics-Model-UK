from modules import input_0D, input_barrages, lagoon_operation
from support.auxillary import *
from support import tools
import numpy as np
import copy


# initialise the inputs for all lagoons in being modelled and store in a dictionary

def custom_generate_inputs(location_list, shift, area=1, operation='two-way'):
    '''
    Function generated for australia case
    :param location_list: point_0 and point_1
    :param run_number:
    :param shift:
    :return:
    '''


    inputs_dict = {}
    input_directory = "inputs/"

    for location in location_list:
        print(location)
        inputs_dict[location] = {}
        # elevation_time_series = input_0D.read_outer_elevations('inputs_0D/WL_' + run_number + '_' + location + ".npy")
        elevation_time_series = input_0D.read_outer_elevations(input_directory + location + ".npy")
        area_elevation_curve = input_0D.lagoon_idealised_area_case(area)

        """
            Simulation setup
            """
        # Time step
        Dt = 120.
        # Number of M_2 tidal cycles
        cycles = 704
        # Conversion of cycles in seconds
        t_sim = cycles * 12.42 * 3600

        # initialise control / design parameteres - reading LagoonSpecs.dat
        # Preliminary calculation of energy and capacity.
        print(location)
        print("Assumed area = ", area, "km2")

        tidal_range, potential_energy = determine_mean_tidal_range_and_energy(np.arange(0, t_sim, Dt),
                                                                              elevation_time_series(np.arange(0, t_sim, Dt)),
                                                                              area=area * 1e6)
        print("Mean tidal range (m) =  " + str(tidal_range) + "\n")
        print("Potential Energy (MWh) =" + str(potential_energy) + "\n")
        capacity = determine_capacity(area, tidal_range, efficiency=0.4 ,capacity_factor=0.2)

        """
        Configuration
        """
        H_m, P_c = determine_mean_head(np.arange(0, t_sim, Dt), elevation_time_series(np.arange(0, t_sim, Dt)))
        print("Mean Range = ", H_m, "m")

        # Suggested turbine parameters
        turbine_params = {"f_g": 50, "g_p": 96, "g": 9.807, "t_d": 7.35,
                          "t_cap": P_c, "h_cap": 0.75 * H_m, "dens": 1025, "h_min": 1.00,
                          "eta": [0.93, 0.83], "options": 1}

        H_m, P_c = determine_mean_head(np.arange(0, t_sim, Dt), elevation_time_series(np.arange(0, t_sim, Dt)),
                                       turbine_params=turbine_params)

        print("Rated Head = ", 0.75 * H_m, "m, Turbine Capacity = ", P_c, " MW")

        turbines = capacity / P_c
        sluices = turbines / 2

        original_control, lagoon_params = input_barrages.input_predefined_barrage_specs(turbines, sluices,
                                                                                        operation=operation,
                                                                                        turbine_parameters=turbine_params)

        lagoon_status = input_barrages.initialise_barrage(len(original_control))
        status, params = {}, {}
        status.update(lagoon_status[0])
        params.update(lagoon_params[0])
        Dt = 100
        cycles_ramp = 10
        simulation_ramp = {"t": (cycles_ramp * 12.42 * 3600) + shift, "Dt": Dt, "start_t": - cycles_ramp * 12.42 * 3600}
        status["m_t"] = - cycles_ramp * 12.42 * 3600
        starting_status = lagoon_operation.tidal_lagoon_0d_model(simulation_ramp, elevation_time_series,
                                                                 area_elevation_curve, status, original_control[0], params)
        inputs_dict[location]['elev_ts'] = elevation_time_series
        inputs_dict[location]['aec'] = area_elevation_curve
        inputs_dict[location]['control_orig'] = original_control[0]
        inputs_dict[location]['status'] = starting_status
        inputs_dict[location]['lag_params'] = lagoon_params

        # print('generating inputs for:  ' + location)
        # print('original control:  ', original_control[0])
        # print('starting status:  ', starting_status)

    return inputs_dict



def generate_inputs(location_list, run_number, shift):

    inputs_dict = {}

    for location in location_list:

        inputs_dict[location] = {}
        elevation_time_series = input_0D.read_outer_elevations('inputs_0D/WL_' + run_number + '_' + location + ".npy")
        area_elevation_curve = input_0D.read_area_elevation_curve('inputs_0D/area_' + location + ".npy",
                                                                  depth_correction= 0.)
        # original_control, lagoon_params = input_barrages.input_barrage('inputs_0D/LagoonSpecs_' + location + ".dat")
        original_control, lagoon_params = input_barrages.input_barrage('inputs_0D/LagoonSpecs_' + location + "_2hrs.dat")
        lagoon_status = input_barrages.initialise_barrage()
        status, params = {}, {}
        status.update(lagoon_status[0])
        params.update(lagoon_params[0])
        Dt = 100
        cycles_ramp = 10
        simulation_ramp = {"t": (cycles_ramp * 12.42 * 3600) + shift, "Dt": Dt, "start_t": - cycles_ramp * 12.42 * 3600}
        starting_status = lagoon_operation.tidal_lagoon_0d_model(simulation_ramp, elevation_time_series,
                                                                 area_elevation_curve, status, original_control[0], params)
        inputs_dict[location]['elev_ts'] = elevation_time_series
        inputs_dict[location]['aec'] = area_elevation_curve
        inputs_dict[location]['control_orig'] = original_control[0]
        inputs_dict[location]['status'] = starting_status
        inputs_dict[location]['lag_params'] = lagoon_params

        # print('generating inputs for:  ' + location)
        # print('original control:  ', original_control[0])
        # print('starting status:  ', starting_status)

    return inputs_dict


# the objective function for income - all lagoons optimised together
def obj_func_baseload(x, simulation_main, input_dict, i, lagoon_list, shift, obj_func='minimum_energy'):
    power_sum = np.zeros(int(simulation_main['t'] / simulation_main['Dt']))
    total_energy = 0
    for s, lagoon in enumerate(lagoon_list):
        ctrl = copy.deepcopy(input_dict[lagoon]['control_orig'])
        status_copy = input_dict[lagoon]['status'].copy()
        z = [0] * i
        ctrl["h_t"][0], ctrl["h_t"][1] = x[2*s], x[2*s+1]
        z.insert(i, ctrl)
        # print(lagoon, ctrl)
        lag_stat, output = lagoon_operation.tidal_lagoon_0d_model(simulation_main, input_dict[lagoon]['elev_ts'],
                                                              input_dict[lagoon]['aec'], status_copy, z,
                                                              input_dict[lagoon]['lag_params'][0],
                                                              export_output=True, adaptive_operation=z, shift=shift)
        power_gen = output[:, 4].copy()
        power_sum = power_sum + power_gen
        total_energy += np.sum(power_gen * simulation_main["Dt"] / 3600)
    minimum_power = np.min(power_sum)
    readings_over_0 = np.sum(power_sum > 0)  # want to be higher
    seconds_over_0 = readings_over_0*simulation_main['Dt']
    if obj_func == 'minimum_power':
        print(str(minimum_power) + 'MW')
        return -minimum_power
    elif obj_func == 'generation_time':
        print(str(seconds_over_0) + ' seconds')
        return -seconds_over_0
    else:
        TypeError('Enter valid objective function')


# only used when optimising individual schemes to maximise energy
def obj_func_energy_multiple(x, simulation_main, input_dict, i, lag_list, shift):
    total_energy = 0
    for s, lagoon in enumerate(lag_list):
        ctrl = copy.deepcopy(input_dict[lagoon]['control_orig'])
        status_copy = input_dict[lagoon]['status'].copy()
        z = [0] * i
        simulation_main["t"] = 2 * 12.42 * 3600
        ctrl["h_t"][0], ctrl["h_t"][1] = x[2 * s], x[2 * s + 1]
        z.insert(i, ctrl)
        z.append(input_dict[lagoon]['control_orig'])
        lagoon_stat = lagoon_operation.tidal_lagoon_0d_model(simulation_main, input_dict[lagoon]['elev_ts'],
                                                             input_dict[lagoon]['aec'], status_copy, z,
                                                             input_dict[lagoon]['lag_params'][0],
                                                             adaptive_operation=z, shift=shift)
        total_energy += lagoon_stat['E']
    simulation_main["t"] = 12.42 * 3600
    print(str(total_energy) + 'Mwh per double tidal cycle')
    return -total_energy


# the objective function for optimising energy - one lagoon passed through at a time
def obj_func_energy_single(x, simulation_main, lagoon_dict, i, shift, operation):
    ctrl = copy.deepcopy(lagoon_dict['control_orig'])
    status_copy = lagoon_dict['status'].copy()
    z = [0] * i
    simulation_main["t"] = 2 * 12.42 * 3600

    if operation == 'ebb':
       ctrl["h_t"][0] = x[0]
    elif operation == 'two-way':
        ctrl["h_t"][0], ctrl["h_t"][1] = x[0], x[1]
    elif operation == 'two-way-pump':
        ctrl["h_t"][0], ctrl["h_t"][1], ctrl["t_p"][0], ctrl["t_p"][1] = x[0], x[1], x[2], x[3]

    z.insert(i, ctrl)
    z.append(lagoon_dict['control_orig'])
    lag_stat = lagoon_operation.tidal_lagoon_0d_model(simulation_main, lagoon_dict['elev_ts'], lagoon_dict['aec'],
                                                      status_copy, z, lagoon_dict['lag_params'][0],
                                                      adaptive_operation=z, shift=shift)
    energy = lag_stat["E"]
    simulation_main["t"] = 12.42 * 3600
    print(ctrl, energy)
    return -energy


def gauge_func_status(x, simulation_main, lagoon_dict, operation='two-way'):
    ctrl = copy.deepcopy(lagoon_dict['control_orig'])
    status_copy = lagoon_dict['status'].copy()

    if operation == 'ebb':
        ctrl["h_t"][0] = x[0]
    elif operation == 'two-way':
        ctrl["h_t"][0], ctrl["h_t"][1] = x[0], x[1]
    elif operation == 'two-way-pump':
        ctrl["h_t"][0], ctrl["h_t"][1], \
        ctrl["t_p"][0], ctrl["t_p"][1] = x[0], x[1], x[2], x[3]

    lagoon_status = lagoon_operation.tidal_lagoon_0d_model(simulation_main, lagoon_dict['elev_ts'], lagoon_dict['aec'],
                                                           status_copy, ctrl, lagoon_dict['lag_params'][0])
    return lagoon_status


def gauge_func_baseload(x_res, simulation_main, input_dict, lagoon_list, shift, to_return):
    power_sum = np.zeros(int(simulation_main['t'] / simulation_main['Dt']))
    total_energy = 0
    for s, lagoon in enumerate(lagoon_list):
        print(lagoon)
        ctrl = copy.deepcopy(input_dict[lagoon]['control_orig'])
        status_copy = input_dict[lagoon]['status'].copy()
        ctrl["h_t"][0], ctrl["h_t"][1] = x_res[2 * s], x_res[2 * s + 1]
        lag_stat, output = lagoon_operation.tidal_lagoon_0d_model(simulation_main, input_dict[lagoon]['elev_ts'],
                                                                  input_dict[lagoon]['aec'], status_copy, ctrl,
                                                                  input_dict[lagoon]['lag_params'][0],
                                                                  export_output=True, shift=shift)
        power_gen = output[:, 4].copy()
        power_sum = power_sum + power_gen
        total_energy += np.sum(power_gen * simulation_main["Dt"] / 3600)
    minimum_power = np.min(power_sum)
    minimum_energy = minimum_power * simulation_main["Dt"] / 3600
    readings_over_0 = np.sum(power_sum > 0)  # want to be higher
    total_readings = power_sum.shape[0]
    if to_return == 'minimum_power':
        return minimum_power
    elif to_return == 'generation_time':
        return int(readings_over_0), int(total_readings)
    elif to_return == 'total_energy':
        return total_energy
    else:
        TypeError('Enter output option')
