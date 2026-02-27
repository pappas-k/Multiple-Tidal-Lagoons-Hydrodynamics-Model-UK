import uptide
from scipy.interpolate import interp1d
import numpy as np
import sys

sys.path.append('../')

from tools.peaks import peakdet
import datetime


def signal_reconstruction(time, series, constituents, signal_time, t_start=datetime.datetime(2003, 5, 6, 8, 0), dt=0.0):
    """
    Conducts signal reconstruction using an array of time and a given series, the selected tide constituents
    and outputs a signal for a desired fixed period.

    :param time: time array
    :param series: series array (can be elevation or velocity component)
    :param constituents: constituents to consider in the harmonic analysis
    :param signal_time: array for the time to reconstruct (in seconds)
    :param t_start: starting time of reconstructed signal (datetime object)
    :param dt: difference in seconds from start time
    :return:
    """
    tide = uptide.Tides(constituents)
    tide.set_initial_time(t_start)
    amplitude, phase = uptide.harmonic_analysis(tide, series, time)
    return interp1d(signal_time, tide.from_amplitude_phase(amplitude, phase, signal_time + dt))


def produce_reconstructed_signal(time, series, constituents, t_start=datetime.datetime(2003, 5, 6, 8, 0),
                                 start_time=0.0,
                                 end_time=3600 * 24 * 30, increment=100):
    """
    Produces reconstructed signal based on a time series and a series of observed data (elevation, velocity components)

    :param time: time array
    :param series: series array (can be elevation or velocity component)
    :param constituents: constituents to consider in the harmonic analysis
    :param t_start: starting time of reconstructed signal (datetime object)
    :param start_time: start time of reconstructed signal
    :param end_time: end time of reconstructed signal
    :param increment: increment in reconstructed signal
    :return:
    """
    n_increment = int((end_time - start_time) / increment + 1)
    t_0d = np.linspace(start_time, end_time, int(n_increment))
    index = np.arange(n_increment)
    signal = signal_reconstruction(time, series, constituents, t_0d, t_start)
    signal_predicted = np.zeros(n_increment)
    signal_predicted[:] = signal(start_time + index[:] * increment)
    return t_0d, signal_predicted


def determine_mean_tidal_range(t, eta, area=1e6):
    """
    Function that quantifies the range and energy from an elevation signal

    :param t: time array
    :param eta: elevation array
    :param area: area occupied for potential energy calculation
    :return: average range, potential energy.
    """
    dens, grav = 1025, 9.81
    Dt = 0.1

    """
    eta will be our water level series
    t will be our time series
    """

    t = t.tolist()
    eta = eta.tolist()
    H, t_p = [], []

    """
    In the following lines we extract the water elevation peaks (HW,LW)
    """
    maxtab, mintab = peakdet(eta, Dt, t)
    for i in range(max(len(maxtab), len(mintab))):
        try:
            t_p.append(np.average([maxtab[i, 0], mintab[i, 0]]))
            H.append(abs(maxtab[i, 1] - mintab[i, 1]))
        except(IndexError):
            pass
        try:
            t_p.append(np.average([maxtab[i + 1, 0], mintab[i, 0]]))
            H.append(abs(maxtab[i + 1, 1] - mintab[i, 1]))
        except(IndexError):
            pass

    H = np.array(H)
    E = np.zeros(len(H))
    E[:] = 0.5 * grav * dens * area * H[:] * H[:] / (1e3 * 3600)

    """
    Here we calculate the Energy for each transition from HW to LW (converted to kWh)
    """
    range1 = np.average(H)
    # print(range1)
    # energy = 1411 * 0.5 * grav * dens * area * range1 * range1 /(1e3 * 3600)
    energy = np.sum(E) * 12
    return range1, energy


"""
TIDAL STREAM RELATED FUNCTIONS
"""


def theoretical_tidal_stream_power(u, density=1025, diameter=20):
    """
    :param u: inflow velocity (m/s)
    :param density: medium density (kg/m^3)
    :return: power output (MW)
    """
    return 1 / 2 * (np.pi * diameter ** 2) / 4 * density * u ** 3


# Thrust coefficient Martin-Short (2015)
def c_t(u, u_in=1.0, u_r=2.5, c_t=0.6):
    """
    :param u: velocity (m/s)
    :param u_in: cut-in speed (default:1.0 m/s)
    :param u_r: rated velocity (m/s)
    :return: c_t
    """
    return np.where(u <= u_in, 0, np.where(u <= u_r, c_t, c_t * u_r ** 3 / u ** 3))


# power coefficient calculation Martin-Short (2015)
def c_p(u):
    """
    :param u: velocity (m/s)
    :return: power coefficient
    """
    return 1 / 2 * (1 + np.sqrt(1 - c_t(u))) * c_t(u)


# Turbine power calculation
def p_t(u, diameter=20, density=1025, capacity=1.5e6):
    """
    :param u: inflow velocity (m/s)
    :param diameter: turbine diameter (m)
    :param density: medium density (kg/m^3)
    :return: power output (MW)
    """
    p = 1 / 2 * density * c_p(u) * (np.pi * (diameter / 2) ** 2 * u ** 3)
    return np.where(p <= capacity, p, capacity)


# Force calculation
def f_t(u, diameter=20, density=1025, *kwargs):
    """
    :param u: velocity magnitude (m/s)
    :param diameter: turbine diameter (m)
    :param density: fluid density (kg/m^3)
    :param kwargs: other keyword args
    :return:
    """
    return 1 / 2 * density * c_t(u) * np.pi * (diameter / 2) ** 2 * u ** 2
