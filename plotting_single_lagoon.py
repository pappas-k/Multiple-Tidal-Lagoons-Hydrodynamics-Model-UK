# Plot the water elevations (outer and inner of barrage) and the power output
# of a lagoon in the list ['SW','CA','WA','CO','LI','BL','SO']
# ================================================================================
import h5py
import matplotlib.pyplot as plt
import numpy as np


lag = 'SW'

df = h5py.File(f'outputs/outputs_run/diagnostic_lagoon_{lag}.hdf5', 'r')
print(df.items())

# Index library for output
index_conv = {"t": 0, "h_o": 1, "h_i": 2, "DZ": 3, "P": 4, "E": 5,
              "m": 6, "Q_t": 7, "Q_s": 8, "m_dt": 9, "m_t": 10, "f_r": 11}

fig, ax = plt.subplots(2, figsize=(10, 4), sharex="all")

# Plotting Elevations in time (black for outer water levels and red for inner)
ax[0].plot(df["operation_output"][:, index_conv["t"]], df["operation_output"][:, index_conv["h_o"]],
           color="black", label='$\eta_o$')
ax[0].plot(df["operation_output"][:, index_conv["t"]], df["operation_output"][:, index_conv["h_i"]], color="r", label='$\eta_i$ (2D)')
ax[0].set_ylabel("$\\eta$ (m)")

# Plotting Power output in time
ax[1].plot(df["operation_output"][:, index_conv["t"]], df["operation_output"][:, index_conv["P"]], color="red")
ax[1].set_ylabel("$P$ (MW)")
ax[1].set_xlabel("$t$ (sec)")
ax[0].legend()
fig.subplots_adjust(hspace=0)
plt.show()


