# %% this reads the data and plots it
import h5py as h5
import numpy as np
import matplotlib.pyplot as plt

with h5.File('orbital/Au.h5', 'r') as f:
    indices = list(f.keys())
    indices.sort(key=lambda x: int(x))
    for index in indices:
        r = f[index]['r'][:]
        density = f[index]['density'][:]
        orbital = f[index].attrs['orbital']
        plt.plot(r, density, label=orbital)
plt.xscale('log')
plt.xlabel('r (Bohr)')
plt.ylabel('Radial density (e/Bohr)')
plt.legend()

# %%
