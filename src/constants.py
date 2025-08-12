"""adpated from ase.units

Physical constants and units derived from CODATA for converting
to and from ase internal units.

adpated from ase chemical_symbols and atomic_numbers
"""

from math import pi, sqrt
from typing import Dict

CODATA = {
    # the "original" CODATA version ase used ever since
    # Constants from Konrad Hinsen's PhysicalQuantities module (1986 CODATA)
    # Add the constant pi used to define the mu0 and hbar here for reference
    # as well
    "1986": {
        "_c": 299792458.0,  # speed of light, m/s
        "_mu0": 4.0e-7 * pi,  # permeability of vacuum
        "_Grav": 6.67259e-11,  # gravitational constant
        "_hplanck": 6.6260755e-34,  # Planck constant, J s
        "_e": 1.60217733e-19,  # elementary charge
        "_me": 9.1093897e-31,  # electron mass
        "_mp": 1.6726231e-27,  # proton mass
        "_Nav": 6.0221367e23,  # Avogadro number
        "_k": 1.380658e-23,  # Boltzmann constant, J/K
        "_amu": 1.6605402e-27,
    },  # atomic mass unit, kg
    # CODATA 1998 taken from
    # https://doi.org/10.1103/RevModPhys.72.351
    "1998": {
        "_c": 299792458.0,
        "_mu0": 4.0e-7 * pi,
        "_Grav": 6.673e-11,
        "_hplanck": 6.62606876e-34,
        "_e": 1.602176462e-19,
        "_me": 9.10938188e-31,
        "_mp": 1.67262158e-27,
        "_Nav": 6.02214199e23,
        "_k": 1.3806503e-23,
        "_amu": 1.66053873e-27,
    },
    # CODATA 2002 taken from
    # https://doi.org/10.1103/RevModPhys.77.1
    "2002": {
        "_c": 299792458.0,
        "_mu0": 4.0e-7 * pi,
        "_Grav": 6.6742e-11,
        "_hplanck": 6.6260693e-34,
        "_e": 1.60217653e-19,
        "_me": 9.1093826e-31,
        "_mp": 1.67262171e-27,
        "_Nav": 6.0221415e23,
        "_k": 1.3806505e-23,
        "_amu": 1.66053886e-27,
    },
    # CODATA 2006 taken from
    # https://doi.org/10.1103/RevModPhys.80.633
    "2006": {
        "_c": 299792458.0,
        "_mu0": 4.0e-7 * pi,
        "_Grav": 6.67428e-11,
        "_hplanck": 6.62606896e-34,
        "_e": 1.602176487e-19,
        "_me": 9.10938215e-31,
        "_mp": 1.672621637e-27,
        "_Nav": 6.02214179e23,
        "_k": 1.3806504e-23,
        "_amu": 1.660538782e-27,
    },
    # CODATA 2010 taken from
    # https://doi.org/10.1103/RevModPhys.84.1527
    "2010": {
        "_c": 299792458.0,
        "_mu0": 4.0e-7 * pi,
        "_Grav": 6.67384e-11,
        "_hplanck": 6.62606957e-34,
        "_e": 1.602176565e-19,
        "_me": 9.10938291e-31,
        "_mp": 1.672621777e-27,
        "_Nav": 6.02214129e23,
        "_k": 1.3806488e-23,
        "_amu": 1.660538921e-27,
    },
    # CODATA 2014 taken from
    # http://arxiv.org/pdf/1507.07956.pdf
    "2014": {
        "_c": 299792458.0,
        "_mu0": 4.0e-7 * pi,
        "_Grav": 6.67408e-11,
        "_hplanck": 6.626070040e-34,
        "_e": 1.6021766208e-19,
        "_me": 9.10938356e-31,
        "_mp": 1.672621898e-27,
        "_Nav": 6.022140857e23,
        "_k": 1.38064852e-23,
        "_amu": 1.660539040e-27,
    },
    # CODATA 2018 taken from
    # https://physics.nist.gov/cuu/Constants/index.html
    "2018": {
        "_c": 299792458.0,  # Exact
        "_mu0": 4.0e-7 * pi,  # Exact
        "_Grav": 6.67430e-11,  # +/- 0.000_15e-11
        "_hplanck": 6.62607015e-34,  # Exact
        "_e": 1.602176634e-19,  # Exact
        "_me": 9.1093837015e-31,  # +/- 0.000_000_0028e-31
        "_mp": 1.67262192369e-27,  # +/- 0.000_000_000_51e-27
        "_Nav": 6.02214076e23,  # Exact
        "_k": 1.380649e-23,  # Exact
        "_amu": 1.66053906660e-27,
    },  # +/- 0.000_000_000_50e-27
}


class Units:
    """
    this is the hard-coded CODATA values
    all other units are dynamically derived from these values upon import of the module
    """

    def __init__(self, codata_version: str = "2014") -> None:
        self._units = None
        self._codata_version = codata_version

    def __getattr__(self, name):
        if self._units is None:
            self._units = self.create_units()
        return self._units.get(name)

    def create_units(self) -> Dict[str, float]:
        """
        Creates a dictionary containing physical constants and units.

        Parameters:
        codata_version (str): CODATA version to use.

        Returns:
        Units: Dictionary containing physical constants and units.

        Raises:
        ValueError: If the CODATA version is not implemented.
        """
        try:
            units = CODATA[self._codata_version]
        except KeyError as exc:
            raise ValueError(
                f'CODATA version "{self._codata_version}" not implemented'
            ) from exc
        # derived from the CODATA values
        units["_eps0"] = 1 / units["_mu0"] / units["_c"] ** 2  # permittivity of vacuum
        units["_hbar"] = units["_hplanck"] / (2 * pi)  # Planck constant / 2pi, J s

        units["Ang"] = units["Angstrom"] = 1.0
        units["nm"] = 10.0
        units["Bohr"] = (
            4e10
            * pi
            * units["_eps0"]
            * units["_hbar"] ** 2
            / units["_me"]
            / units["_e"] ** 2
        )  # Bohr radius

        units["eV"] = 1.0
        units["Hartree"] = (
            units["_me"]
            * units["_e"] ** 3
            / 16
            / pi**2
            / units["_eps0"] ** 2
            / units["_hbar"] ** 2
        )
        units["kJ"] = 1000.0 / units["_e"]
        units["kcal"] = 4.184 * units["kJ"]
        units["mol"] = units["_Nav"]
        units["Rydberg"] = 0.5 * units["Hartree"]
        units["Ry"] = units["Rydberg"]
        units["Ha"] = units["Hartree"]

        units["second"] = 1e10 * sqrt(units["_e"] / units["_amu"])
        units["fs"] = 1e-15 * units["second"]

        units["kB"] = units["_k"] / units["_e"]  # Boltzmann constant, eV/K

        units["Pascal"] = (1 / units["_e"]) / 1e30  # J/m^3
        units["GPa"] = 1e9 * units["Pascal"]
        units["bar"] = 1e5 * units["Pascal"]

        units["Debye"] = 1.0 / 1e11 / units["_e"] / units["_c"]
        units["alpha"] = (
            units["_e"] ** 2 / (4 * pi * units["_eps0"]) / units["_hbar"] / units["_c"]
        )  # fine structure constant
        units["invcm"] = (
            100 * units["_c"] * units["_hplanck"] / units["_e"]
        )  # cm^-1 energy unit

        # Derived atomic units that have no assigned name:
        units["_aut"] = units["_hbar"] / (
            units["alpha"] ** 2 * units["_me"] * units["_c"] ** 2
        )  # atomic unit of time, s
        units["_auv"] = (
            units["_e"] ** 2 / units["_hbar"] / (4 * pi * units["_eps0"])
        )  # atomic unit of velocity, m/s
        units["_auf"] = (
            units["alpha"] ** 3 * units["_me"] ** 2 * units["_c"] ** 3 / units["_hbar"]
        )  # atomic unit of force, N
        units["_aup"] = (
            units["alpha"] ** 5
            * units["_me"] ** 4
            * units["_c"] ** 5
            / units["_hbar"] ** 3
        )  # atomic unit of pressure, Pa

        units["AUT"] = units["second"] * units["_aut"]  # atomic unit of time

        # SI units
        units["m"] = 1e10 * units["Ang"]  # metre
        units["kg"] = 1.0 / units["_amu"]  # kilogram
        units["s"] = units["second"]  # second
        units["A"] = 1.0 / units["_e"] / units["s"]  # ampere
        units["J"] = units["kJ"] / 1000  # Joule = kg * m**2 / s**2
        units["C"] = 1.0 / units["_e"]  # Coulomb = A * s

        return units


units = Units("2018")

chemical_symbols = [
    # 0
    "X",
    # 1
    "H",
    "He",
    # 2
    "Li",
    "Be",
    "B",
    "C",
    "N",
    "O",
    "F",
    "Ne",
    # 3
    "Na",
    "Mg",
    "Al",
    "Si",
    "P",
    "S",
    "Cl",
    "Ar",
    # 4
    "K",
    "Ca",
    "Sc",
    "Ti",
    "V",
    "Cr",
    "Mn",
    "Fe",
    "Co",
    "Ni",
    "Cu",
    "Zn",
    "Ga",
    "Ge",
    "As",
    "Se",
    "Br",
    "Kr",
    # 5
    "Rb",
    "Sr",
    "Y",
    "Zr",
    "Nb",
    "Mo",
    "Tc",
    "Ru",
    "Rh",
    "Pd",
    "Ag",
    "Cd",
    "In",
    "Sn",
    "Sb",
    "Te",
    "I",
    "Xe",
    # 6
    "Cs",
    "Ba",
    "La",
    "Ce",
    "Pr",
    "Nd",
    "Pm",
    "Sm",
    "Eu",
    "Gd",
    "Tb",
    "Dy",
    "Ho",
    "Er",
    "Tm",
    "Yb",
    "Lu",
    "Hf",
    "Ta",
    "W",
    "Re",
    "Os",
    "Ir",
    "Pt",
    "Au",
    "Hg",
    "Tl",
    "Pb",
    "Bi",
    "Po",
    "At",
    "Rn",
    # 7
    "Fr",
    "Ra",
    "Ac",
    "Th",
    "Pa",
    "U",
    "Np",
    "Pu",
    "Am",
    "Cm",
    "Bk",
    "Cf",
    "Es",
    "Fm",
    "Md",
    "No",
    "Lr",
    "Rf",
    "Db",
    "Sg",
    "Bh",
    "Hs",
    "Mt",
    "Ds",
    "Rg",
    "Cn",
    "Nh",
    "Fl",
    "Mc",
    "Lv",
    "Ts",
    "Og",
]

atomic_numbers = {}
for Z, symbol in enumerate(chemical_symbols):
    atomic_numbers[symbol] = Z
