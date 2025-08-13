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
elements = [6,7,14, 79]  # C, N, Si, Au
g_values = np.geomspace(1e-2,1.5,100)
# g_values = np.linspace(0,20,100)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

for Z in elements:
    sf = ScatteringFactor(Z=Z)
    
    # Calculate form factors
    fxg = sf.xray_form_factor(g=g_values, mode='total')
    feg = sf.electron_form_factor(g=g_values, mode='total')
    
    # Plot
    ax1.plot(g_values*2, fxg, '-', linewidth=2, label=f'{sf.element} (Z={Z})')
    ax2.plot(g_values*2, feg, '-', linewidth=2, label=f'{sf.element} (Z={Z})')

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
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from abtem.parametrizations import LobatoParametrization, KirklandParametrization
from ase.io import read

import abtem

abtem.config.set({"diagnostics.progress_bar": False});
sys.path.insert(0, os.path.abspath('/home/zzhang/OneDrive/code/coredensity') )

from src.scattering_factor import ScatteringFactor
symbols = ["C", "Si", 'Cu', "Au",'U']


lobato = LobatoParametrization()

potentials = lobato.line_profiles(symbols, cutoff=2, name="potential")
lobato = lobato.line_profiles(
    symbols, cutoff=3, name="scattering_factor"
)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4))
ax1.grid(True, alpha=0.3)
ax2.grid(True, alpha=0.3)

visualization = potentials.show(ax=ax1, legend=False)
visualization.set_ylim([-1e2, 2e3])

lobato.show(legend=True, ax=ax2);

kirkland = KirklandParametrization()

potentials = kirkland.line_profiles(symbols, cutoff=2, name="potential")
kirkland = kirkland.line_profiles(
    symbols, cutoff=3, name="scattering_factor"
)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4))
ax1.grid(True, alpha=0.3)
ax2.grid(True, alpha=0.3)

visualization = potentials.show(ax=ax1, legend=False)
visualization.set_ylim([-1e2, 2e3])

lobato.show(legend=True, ax=ax2);
"""Compare scattering factors for different elements"""
elements = [6,14, 29, 79,92]  # C, Si, Cu, Au, U
g_values = np.geomspace(1e-2,1.5,100)
# g_values = np.linspace(0,20,100)
k = np.arange(lobato.sampling, 3+lobato.sampling, lobato.sampling)
 
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
for Z in elements:
    sf = ScatteringFactor(Z=Z)
    
    # Calculate form factors
    fxg = sf.xray_form_factor(g=g_values, mode='total')
    feg = sf.electron_form_factor(g=g_values, mode='total')
    
    # Plot
    ax1.plot(g_values*2, fxg, '-', linewidth=2, label=f'Dirac {sf.element} (Z={Z})')
    ax2.plot(g_values*2, feg, '-', linewidth=2, label=f'Dirac {sf.element} (Z={Z})')
    i = elements.index(Z)
    ax2.plot(k, lobato.array[i], '--', label=f'abTEM Lobato {sf.element} (Z={Z})')
    ax2.plot(k, kirkland.array[i], ':', label=f'abTEM Kirkland {sf.element} (Z={Z})')

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
