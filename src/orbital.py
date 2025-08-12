from ast import Tuple
import re
from abc import ABCMeta
from typing import List,Tuple
import numpy as np

from .constants import chemical_symbols
from .utils import config_tuples_to_config_str, orbital_configuration


class Orbital(metaclass=ABCMeta):
    """
    Class for storing information about an orbital.

    Parameters:
        Z (int): The atomic number of the element.
        n (int): The principal quantum number.
        l (int): The angular momentum quantum number.
        spin (float): The spin quantum number.
        core_hole (bool): Whether the orbital is a core hole.

    Attributes:
        Z (int): The atomic number of the element.
        n (int): The principal quantum number.
        l (int): The angular momentum quantum number.
        spin (float): The spin quantum number.
        core_hole (bool): Whether the orbital is a core hole.
        element (str): The element of the subshell.
        edge_name (str): The edge name of the subshell.
        orbital_name (str): The orbital name of the subshell.
        config_str (str): The electronic configuration of the subshell.
        electron_configuration (list of tuples): The electronic configuration of the subshell.
        subshell_index (int): The orbital index of the subshell.
        subshell_occupancy (int): The number of electron in the subshell.
        spin_occupancy_ratio (float): The spin-orbital occupancy ratio of the subshell.
    """

    def __init__(
        self,
        Z: int = 1,
        n: int = 1,
        l: int = 0,
        spin: float = 0.5,
        core_hole: bool = False,
        num_unoccupied_states: int = 10,
    ) -> None:
        self._Z = Z
        self._n = n
        self._l = l
        self._spin = spin
        self.core_hole = core_hole
        self._config = None
        self.num_unoccupied_states = num_unoccupied_states
        self._unoccupied_states = None

    @property
    def Z(self):
        """
        Returns the atomic number of the element.
        """
        return self._Z

    @property
    def n(self):
        """
        Returns the principal quantum number.

        Returns:
            int: the principal quantum number.
        """
        return self._n

    @property
    def l(self):
        """
        Returns the angular momentum quantum number.

        Returns:
            int: the angular momentum quantum number.
        """

        return self._l

    @property
    def spin(self):
        """spin quantum number

        Returns:
            float: spin quantum number
        """
        if self.l == 0:
            self._spin = 1 / 2
        return self._spin

    @property
    def j(self):
        """total angular momentum quantum number

        Returns:
            float: total angular momentum quantum number
        """
        return self.l + self.spin

    @property
    def kappa(self):
        """relativistic quantum number

        Returns:
            int: relativistic quantum number
        """
        return np.array(-self.spin * (2 * self.j + 1)).astype(int)

    @property
    def subshell_index(self):
        """
        Returns the orbital index of the subshell.

        Returns:
            int: the orbital index of the subshell.
        """
        config_tuples = self.electron_configuration
        subshell_index = [shell[:2] for shell in config_tuples].index((self.n, self.l))
        return subshell_index

    @property
    def subshell_occupancy(self):
        """
        Returns the number of electron in the subshell.

        Returns:
            int: the number of electron in the subshell.
        """
        config_tuples = self.electron_configuration
        subshell_occupancy = config_tuples[self.subshell_index][-1]
        return subshell_occupancy

    @property
    def spin_occupancy_ratio(self):
        """
        Returns the spin-orbital occupancy ratio of the subshell.
        Defined as the number of target spin electrons divided by the total number of electrons in the current subshell (can be open or closed).

        Returns:
            float: the spin-orbital occupancy ratio of the subshell.
        """
        # return np.array((2 * self.j + 1) / (4 * self.l + 2))
        full_shell = int(4 * self.l + 2)
        occ = self.subshell_occupancy
        if self.subshell_occupancy == full_shell:
            return np.array((2 * self.j + 1) / (4 * self.l + 2))
        else:
            spin_down = min(occ, 2 * self.l)
            spin_up = max(0, occ - 2 * self.l)
            spin_selected = spin_down if self.spin == -1 / 2 else spin_up
            return spin_selected / occ

    @property
    def element(self):
        """
        Returns the element of the subshell.

        Returns:
            str: the element of the subshell.
        """
        return chemical_symbols[self.Z]

    @property
    def edge_name(self):
        """
        Returns the edge name of the subshell.

        Returns:
            str: the edge name of the subshell.
        """
        assert self.n > 0, "The edge name is only defined for the bound state"
        edge = ["K", "L", "M", "N", "O", "P", "Q"][self.n - 1]
        if self.spin == -1 / 2:
            idx = self.l * 2
        else:
            idx = self.l * 2 + 1
        # if self.n == 1:
        #     return edge
        # else:
        return edge + str(idx)

    @property
    def orbital_name(self):
        """
        Returns the orbital name of the subshell.

        Returns:
            str: the orbital name of the subshell.
        """
        orbital = [
            "s",
            "p",
            "d",
            "f",
            "g",
            "h",
            "i",
            "k",
            "l",
            "m",
            "n",
            "o",
            "q",
            "r",
            "s",
            "t",
            "u",
            "v",
            "w",
            "x",
            "y",
            "z",
        ][self.l]
        if self.spin == -1 / 2:
            spin = "down"
        else:
            spin = "up"
        return str(self.n) + orbital + "_" + spin

    @property
    def orbital_symbol(self):
        """
        Returns the orbital name of the subshell.

        Returns:
            str: the orbital name of the subshell.
        """
        orbital = [
            "s",
            "p",
            "d",
            "f",
            "g",
            "h",
            "i",
            "k",
            "l",
            "m",
            "n",
            "o",
            "q",
            "r",
            "s",
            "t",
            "u",
            "v",
            "w",
            "x",
            "y",
            "z",
        ][self.l]
        spin = f"{int(2*self.j)}/2"
        return str(self.n) + orbital + spin

    @property
    def config_str(self):
        """electronic configuration in string

        Returns:
            str: electronic configuration in string
        """
        config_tuples = self.electron_configuration
        config = config_tuples_to_config_str(config_tuples)
        return config

    @property
    def electron_configuration(self):
        """
        Returns the electronic configuration of the subshell.

        Returns:
            list of tuples: the electronic configuration of the subshell.
        """
        if self._config is None:
            self._config = self.calculate_electron_configuration()
        return self._config

    def set_electron_configuration(self, config_tuples):
        """
        set the electronic configuration of the subshell.

        Args:
            config_tuples (list of tuples): the electronic configuration of the subshell.
        """
        self._config = config_tuples

    def calculate_electron_configuration(self):
        """
        Calculates the electronic configuration of the subshell based on the orbital configuration and the oxidation state.

        Returns:
            list of tuples: the electronic configuration of the subshell.
        """
        config_tuples = self.orbital_configuration
        return config_tuples

    @property
    def valence_shell(self):
        """
        Returns the valence shell of the atom.

        Returns:
            tuple: the valence shell of the atom.
        """
        valence_shell = self.electron_configuration[-1]
        return valence_shell

    @property
    def orbital_configuration(self):
        """
        Returns the orbital configuration of the subshell.

        Returns:
            list of tuples: the orbital configuration of the subshell.
        """
        config_tuples = orbital_configuration(self.Z)
        return config_tuples

    @property
    def spin_configuration(self):
        """spin configuration of the atom for getting the term symbol
        Here the spin is in the magnetic representation (half spin up and half spin down), not in the spectroscopic representation with spin ratio = (2j+1)/(4l+2)

        Returns:
            list of tuples: the spin orbital configuration of the atom
        """
        spin_orbital_list = []
        for shell in self.electron_configuration:
            n, l, occ = shell
            for spin in [-1 / 2, 1 / 2]:
                if l == 0 and spin == -1 / 2:
                    continue
                j = l + spin
                full_spin = int(2 * j + 1)
                if spin == -1 / 2:
                    spin_occ = min(occ, full_spin)
                else:
                    spin_occ = max(0, occ - (2 * l))
                spin_orbital_list.append((n, l, spin, spin_occ))
        return spin_orbital_list

    def _update(self, Z, n, l, spin):
        """update the atom with new atomic state

        Args:
            Z (int): atomic number
            n (int): principal quantum number
            l (int): orbital quantum number
            spin (float): spin quantum number
        """
        self._Z = Z
        self._n = n
        self._l = l
        self._spin = spin
