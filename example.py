# %% this reads the data and plots it
import h5py as h5
import numpy as np
import matplotlib.pyplot as plt

with h5.File('database/Au.h5', 'r') as f:
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
# Import from the installed coredensity package
import sys
import os
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from src.scattering_factor import ScatteringFactor

"""Compare scattering factors for different elements"""
elements = [1, 6, 10, 26, 79]  # H, C, Ne, Fe, Au
g_values = np.geomspace(1e-3,20,100)
# g_values = np.linspace(0,20,100)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

for Z in elements:
    sf = ScatteringFactor(Z=Z)
    
    # Calculate form factors
    fxg = sf.xray_form_factor(g=g_values, mode='total')
    feg = sf.electron_form_factor(g=g_values, mode='total')
    
    # Plot
    ax1.plot(g_values, fxg, '-', linewidth=2, label=f'{sf.element} (Z={Z})')
    ax2.plot(g_values, feg, '-', linewidth=2, label=f'{sf.element} (Z={Z})')

ax1.set_xlabel('Scattering vector g (Å⁻¹)')
ax1.set_ylabel('X-ray form factor f_x(g)')
ax1.set_title('X-ray Scattering Factors')
ax1.legend()
ax1.grid(True, alpha=0.3)

ax2.set_xlabel('Scattering vector g (Å⁻¹)')
ax2.set_ylabel('Electron form factor f_e(g) (Å)')
ax2.set_title('Electron Scattering Factors')
ax2.legend()
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('element_comparison.png', dpi=150, bbox_inches='tight')
print("Element comparison saved as 'element_comparison.png'")
plt.show()
    
# %%
