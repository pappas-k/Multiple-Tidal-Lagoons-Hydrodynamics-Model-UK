import h5py
import matplotlib.pyplot as plt
import numpy as np
from inputs.simulation_parameters import *  # af


def plot_data_from_hdf5(filename):
    # Open the HDF5 file and access its contents
    df = h5py.File(filename, 'r')

    # Create a 2x1 subplot grid to display the plots side by side
    fig, ax = plt.subplots(2, figsize=(10, 4), sharex="all")

    # Plot water level data on the first subplot
    ax[0].plot(df["operation_output"][:, index_conv["t"]], df["operation_output"][:, index_conv["h_o"]],
               color="black", label='$\eta_o$')
    ax[0].plot(df["operation_output"][:, index_conv["t"]], df["operation_output"][:, index_conv["h_i"]], color="r", label='$\eta_i$ (2D)')
    ax[0].set_ylabel("$\\eta$ (m)")

    # Plot power output data on the second subplot
    ax[1].plot(df["operation_output"][:, index_conv["t"]], df["operation_output"][:, index_conv["P"]], color="red")
    ax[1].set_ylabel("$P$ (MW)")
    ax[1].set_xlabel("$t$ (sec)")

    # Add a legend to the first subplot
    ax[0].legend()

    # Adjust the spacing between the subplots
    fig.subplots_adjust(hspace=0)

    # Display the plots
    plt.show()

# Index library for output
index_conv = {"t": 0, "h_o": 1, "h_i": 2, "DZ": 3, "P": 4, "E": 5,
              "m": 6, "Q_t": 7, "Q_s": 8, "m_dt": 9, "m_t": 10, "f_r": 11}

# List of filenames
#file_names = ['diagnostic_lagoon_0_-1.hdf5', 'diagnostic_lagoon_1_-1.hdf5', 'diagnostic_lagoon_2_-1.hdf5', ...]
#file_names = [f"outputs/outputs_ramp/diagnostic_lagoon_{i}_-1.hdf5" for i in range(7)]
#file_names = [f"outputs/outputs_run/diagnostic_lagoon_{i}_0.hdf5" for i in range(7)]
#file_names = [f"outputs/outputs_ramp/diagnostic_lagoon_{lag}_-1.hdf5" for lag in lagoons]
file_names = [f"outputs/outputs_run/diagnostic_lagoon_{lag}_0.hdf5" for lag in lagoons]


# Plot data from each file
for filename in file_names:
    plot_data_from_hdf5(filename)