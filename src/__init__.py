"""
Core Density - Electron density and scattering factor calculations

This package provides electron density calculations for all elements,
calculated by solving the Dirac equation, along with X-ray and electron
scattering factor calculations.
"""

from .scattering_factor import ScatteringFactor
from .orbital import Orbital
from .constants import units, chemical_symbols, atomic_numbers

__version__ = "0.0.1"
__author__ = "Zezhong Zhang"

__all__ = [
    'ScatteringFactor',
    'Orbital', 
    'units',
    'chemical_symbols',
    'atomic_numbers'
]