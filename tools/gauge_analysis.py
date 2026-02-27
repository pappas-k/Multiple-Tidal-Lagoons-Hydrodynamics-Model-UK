import sys                                  #af
import os
sys.path.append('../')                      #af
from inputs.simulation_parameters import i_tidegauge_file, run_output_folder, i_constituents, \
    s_year, s_month, s_day, s_hour, s_min, i_dt
import matplotlib
matplotlib.use('Agg')
import h5py
import uptide
import datetime
import matplotlib.pyplot as plt
import numpy as np
from math import sqrt

# Access simulation output file and setup output folder
os.chdir('../')
det_file = run_output_folder + '/diagnostic_detectors_gauges.hdf5'
output_folder = 'data/'
df = h5py.File(det_file, 'r+')

# Constituents available in the data
data_constituents = ['Q1', 'O1', 'P1', 'S1', 'K1', '2N2', 'MU2', 'N2', 'NU2', 'M2', 'L2', 'T2', 'S2', 'K2', 'M4', 'Z0'] # the ones present in the data
thetis_constituents = i_constituents # the ones you want to analyse for

#  Initialise uptide for harmonic analysis
tide = uptide.Tides(thetis_constituents)
tide.set_initial_time(datetime.datetime(s_year,s_month,s_day,s_hour,s_min))  # make sure this is the same as in the tidal forcing of your run

# Read tidegauge file
tidegauge_file = i_tidegauge_file
gauge_names = np.loadtxt(tidegauge_file, skiprows=1, usecols=(0,), dtype=str, delimiter=',')
gauge_amps = np.loadtxt(tidegauge_file, skiprows=1, usecols=np.arange(7, 36, 2),delimiter=',')
gauge_phases = np.loadtxt(tidegauge_file, skiprows=1, usecols=np.arange(8, 36, 2),delimiter=',')

# Extract time array and setup harmonic analysis timestep based on dt
t = df['time']
dt= i_dt

# Ignore spin up time if necessary
# spin = int(8*24*60*60/dt)
spin = 0


def plot_amplitudes(constituent):
    """
    Calculate amplitude for a particular constituent and compare against known data
    :param constituent: constituent value
    :return:
    """
    cno_data = (np.array(data_constituents)==constituent).argmax()
    #print(cno_data)
    cno_thetis = (np.array(thetis_constituents)==constituent).argmax()
    #print(cno_thetis)
    plt.figure()
    maxamp = 0.0
    rmsesum = []
    guagesum = 0
    for name, data in df.items():
        #print(name)
        name=name.split('_')[0]
        if name not in gauge_names:
            #print('Not plotting detector:', name)
            continue
        ind = (gauge_names==name).argmax()
        amp, pha = uptide.harmonic_analysis(tide, data[spin:,0], t[spin:])
        plt.plot(gauge_amps[ind, cno_data], amp[cno_thetis], '.')
        plt.annotate(name, (gauge_amps[ind, cno_data], amp[cno_thetis]), size=4)
        maxamp = max(maxamp, gauge_amps[ind, cno_data], amp[cno_thetis])
        error = gauge_amps[ind, cno_data] - amp[cno_thetis]
        #print(len(gauge_amps[ind, cno_data]), type(gauge_amps[ind, cno_data]))
        #print(len(amp[cno_thetis]), type(amp[cno_thetis]))
        l2error = (gauge_amps[ind, cno_data] - amp[cno_thetis])**2
        #df.create_dataset('error1', data = error)
        #df.create_dataset('l2error', data = l2error)
        rmsesum.append(l2error)
        guagesum = guagesum + gauge_amps[ind, cno_data]
    rmse = sqrt(sum(rmsesum)/len(rmsesum))
    guageavg = guagesum/len(rmsesum)
    nrmse = rmse/guageavg
    print('{}-constituent rmse = %f'.format(constituent) % rmse)
    print('{}-constituent nrmse = %f'.format(constituent) % nrmse)
    plt.plot([0, maxamp], [0, maxamp], 'k')
    plt.xlabel('Gauge amplitude (m)')
    plt.ylabel('Thetis amplitude (m)')
    plt.title('{}-constituent'.format(constituent))
    plt.savefig(output_folder + '{}-amplitude.png'.format(constituent))
    return rmse, nrmse

# rmse = {'M2':'','O1':'','K1':'','M4':'','S2':''}
# nrmse = {'M2':'','O1':'','K1':'','M4':'','S2':''}

# [rmse['M2'], nrmse['M2']]=plot_amplitudes('M2')
# [rmse['O1'], nrmse['O1']]=plot_amplitudes('O1')
# [rmse['K1'], nrmse['K1']]=plot_amplitudes('K1')
# [rmse['M4'], nrmse['M4']]=plot_amplitudes('M4')
# [rmse['S2'], nrmse['S2']]=plot_amplitudes('S2')

# print(rmse)
# print(nrmse)
#hf.close()
#show()

rmse = {'M2':''}
nrmse = {'M2':''}

text_file_rmse = open(output_folder + "/rmse.txt", "w")
text_file_nrmse = open(output_folder+"/nrmse.txt", "w")


for tc in thetis_constituents:
    #rmse.append(tc,)
    [tc_rmse, tc_nrmse] = plot_amplitudes(tc)
    # print(tc)
    # print(type(tc))
    # print(tc_rmse)
    # print(type(str(tc_rmse)))
    # print(len(str(tc_rmse)))
    rmse[tc] = str(tc_rmse)
    nrmse[tc] = str(tc_nrmse)
    #rmse.update({tc,str(tc_rmse)})
    #nrmse.update({tc,str(tc_nrmse)})
    text_file_rmse.write(tc + ' : ' + str(tc_rmse) + '\n')
    text_file_nrmse.write(tc + ' : ' + str(tc_nrmse) + '\n')

