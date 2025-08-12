import os
import logging
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from scipy.integrate import quad
from abc import ABC
import h5py as h5


from .orbital import Orbital
from .constants import units

ORBITAL_DATABASE_FOLDER = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "database"
)

class ScatteringFactor(Orbital):
    def __init__(self, Z: int = 1, n: int = 1, l: int = 0, spin: float = 0.5) -> None:
        self._fxg = None
        self._feg = None
        super().__init__(Z, n, l, spin)
    
    def load_orbital(self, filename=None):
        with h5.File(filename, "r") as hdf:
            density_frame = pd.DataFrame(
                columns=[
                    "n",
                    "l",
                    "spin",
                    "orbital",
                    "speed",
                    "small_component_ratio",
                    "r",
                    "density_large",
                    "density_small",
                    "density",
                    "wave_large",
                    "wave_small",
                    "r_pot",
                    "vr_pot",
                ]
            )
            indices = list(hdf.keys())
            indices.sort(key=lambda x: int(x))
            for key in indices:
                group = hdf[key]
                n = group.attrs["n"]
                l = group.attrs["l"]
                spin = group.attrs["spin"]
                orbital = group.attrs["orbital"]
                ratio_small = group.attrs["small_component_ratio"]
                speed = group.attrs["speed"]
                r = group["r"][:]
                density_large = group["density_large"][:]
                density_small = group["density_small"][:]
                density = group["density"][:]
                wf_large = group["wave_large"][:]
                wf_small = group["wave_small"][:]
                # save the data to the DataFrame
                new_row = pd.DataFrame(
                    {
                        "n": [n],
                        "l": [l],
                        "spin": [spin],
                        "r": [r],
                        "orbital": [orbital],
                        "speed": [speed],
                        "small_component_ratio": [ratio_small],
                        "density_large": [density_large],
                        "density_small": [density_small],
                        "density": [density],
                        "wave_large": [wf_large],
                        "wave_small": [wf_small],
                    }
                )
                density_frame = pd.concat([density_frame, new_row], ignore_index=True)
        return density_frame
    
    def get_density(self, r = None):
        """get the electron density as a function of the radial grid

        Args:
            r (np.array): radial grid

        Returns:
            pd.DataFrame: density
        """
        filename = os.path.join(ORBITAL_DATABASE_FOLDER, f"{self.element}.h5")

        if not os.path.exists(filename):
            logging.info(f"{filename} does not exist, calculating the orbital data")
        else:
            logging.info(f"loading the orbital data from the {filename} hdf5 file")
        density_frame = self.load_orbital(filename)
        if r is not None:
            density_frame = self.resample_density(r, density_frame)
        return density_frame
    
    def resample_density(self, r, density_frame):
        for idx in density_frame.index:
            r_org = density_frame["r"][idx]
            density = density_frame["density"][idx]
            density_large = density_frame["density_large"][idx]
            density_small = density_frame["density_small"][idx]
            density_interp = interp1d(
                r_org, density, kind="cubic", fill_value=0, bounds_error=False
            )
            density_large_interp = interp1d(
                r_org, density_large, kind="cubic", fill_value=0, bounds_error=False
            )
            density_small_interp = interp1d(
                r_org, density_small, kind="cubic", fill_value=0, bounds_error=False
            )
            density_frame.at[idx, "density"] = density_interp(r)
            density_frame.at[idx, "density_large"] = density_large_interp(r)
            density_frame.at[idx, "density_small"] = density_small_interp(r)
            density_frame.at[idx, "r"] = r
        return density_frame

    def xray_form_factor(
        self, g = np.geomspace(1e-7, 20, 1000), mode="total"
    ):
        """get the x-ray form factor

        Args:
            s (np.array, optional): scattering vector in 1/A. Defaults to np.linspace(1e-7, 20, 1000).

        Returns:
            np.array: x-ray atomic form factor
        """
        assert mode in {"total", "orbital"}, "mode must be 'total' or 'orbital'"
        r = np.geomspace(1e-7, 100, 1000) - 1e-7
        orbital_frame = self.get_density(r=r)
        if mode == "orbital":
            fxg_frame = pd.DataFrame(columns=["n", "l", "spin", "orbital", "s", "fxg"])
            for idx in orbital_frame.index:
                r = orbital_frame["r"][idx]
                density = orbital_frame["density"][idx]
                fxg = self._xray_form_factor(r, density, g)
                # save the data to the DataFrame
                new_row = pd.DataFrame(
                    {
                        "n": orbital_frame["n"][idx],
                        "l": orbital_frame["l"][idx],
                        "spin": orbital_frame["spin"][idx],
                        "orbital": orbital_frame["orbital"][idx],
                        "g": [g],
                        "fxg": [fxg],
                    }
                )
                fxg_frame = pd.concat([fxg_frame, new_row], ignore_index=True)
            return fxg_frame
        else:
            density = orbital_frame["density"].sum()
            fxg = self._xray_form_factor(r, density, g)
            return fxg

    def _xray_form_factor(self, r, density, g):
        """calculate the x-ray form factor

        Args:
            density (np.array): electron density
            g (np.array): scattering vector in 1/A

        Returns:
            np.array: x-ray atomic form factor
        """
        if self._fxg is None:
            density_interp = interp1d(r, density, kind="cubic", fill_value=0, bounds_error=False)
            grind = g * np.pi * 2
            fxg = np.zeros_like(g)
            
            # The function to integrate, with g as an additional argument
            def integrand(r, g_val):
                # Handle the case where g_val * r is zero to avoid division by zero
                if g_val * r == 0:
                    return density_interp(r) # The limit of sin(x)/x as x->0 is 1
                return density_interp(r) * np.sin(g_val * r) / (g_val * r)
            
            # Loop through each value in the grind array and perform the integration
            for i, k in enumerate(grind):
                # Perform the integration using quad
                result, error = quad(
                    integrand,
                    a=0,        # Lower limit of integration
                    b=100,      # Upper limit of integration
                    args=(k,),  # Pass g as an argument to the integrand function
                    epsabs=1e-10, # Absolute error tolerance
                    epsrel=1e-10, # Relative error tolerance
                )

                fxg[i] = result
                print(f"For g = {g[i]:.2f}, the result is {result:.4f} with an error of {error:.2e}")

            
            # Ensure form factor doesn't exceed atomic number
            fxg = np.clip(fxg, 0, self.Z)
            if g[0] == 0:
                # For g=0, the form factor is equal to the atomic number
                fxg[0] = self.Z
                # fxg = np.insert(fxg,0,self.Z) # Insert Z at the beginning for g=0
            self._fxg = fxg
        else:
            fxg = self._fxg
        return fxg

    def electron_form_factor(
        self, g = np.linspace(1e-7, 20, 1000), mode="total"
    ):
        """get the electron form factor

        Args:
            s (np.array, optional): scattering vector in 1/A. Defaults to np.linspace(1e-7, 20, 1000).

        Returns:
            np.array: electron atomic form factor
        """
        assert mode in {"total", "orbital"}, "mode must be 'total' or 'orbital'"
        if mode == "orbital":
            df = self.xray_form_factor(g=g, mode="orbital")
            df["feg"] = None
            for idx in df.index:
                fxg = df["fxg"][idx]
                feg = self.mott_bethe_formula(g, fxg)
                df.at[idx, "feg"] = feg
            return df
        else:
            fxg = self.xray_form_factor(g=g, mode="total")
            feg = self.mott_bethe_formula(g, fxg)
            return feg

    def mott_bethe_formula(self, g, fxg):
        """get the electron form factor

        Args:
            g (np.array, optional): scattering vector in 1/A. Defaults to np.linspace(1e-7, 20, 1000).

        Returns:
            np.array: electron atomic form factor
        """
        # use Mott–Bethe formula
        # https://en.wikipedia.org/wiki/Mott%E2%80%93Bethe_formula
        feg = (self.Z - fxg) / g**2 / (8 * np.pi**2 * units.Bohr)
        max_idx = np.argmax(feg)
        feg[:max_idx] = feg[max_idx]
        return feg
