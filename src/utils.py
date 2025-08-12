import os
from .constants import chemical_symbols



def _set_path(path):
    """
    Internal function to set the parametrization data directory.
    """
    _ROOT = os.path.abspath(os.path.dirname(__file__))
    return os.path.join(_ROOT, path)


def load_orbital_configuration():
    """
    load electronic structure configuration from file
    Returns:
        configurations: electronic structure configuration
    """
    configurations = {}
    with open(_set_path("electron_configurations.txt")) as f:
        for i, line in enumerate(f):
            line = line.strip()
            prefix_start = line.find("[")
            prefix_end = line.find("]")
            if prefix_start > -1:
                line = (
                    configurations[line[prefix_start + 1 : prefix_end]]
                    + line[prefix_end + 1 :]
                )

            configurations[chemical_symbols[i + 1]] = line
    return configurations


def config_str_to_config_tuples(config_str):
    """
    convert electronic structure configuration to tuples
    Returns:
        config_tuples: electronic structure configuration tuples
    """
    azimuthal_number = {"s": 0, "p": 1, "d": 2, "f": 3, "g": 4, "h": 5, "i": 6}
    config_tuples = []
    for subshell_string in config_str.split(" "):
        config_tuples.append(
            (
                int(subshell_string[0]),
                azimuthal_number[subshell_string[1]],
                int(subshell_string[2:]),
            )
        )
    return config_tuples


def orbital_configuration(Z):
    config_tuples = config_str_to_config_tuples(
        load_orbital_configuration()[chemical_symbols[Z]]
    )
    return config_tuples


def config_tuples_to_config_str(config_tuples):
    """
    convert electronic structure configuration to string
    Returns:
        config_str: electronic structure configuration string
    """
    config_str = [
        "".join(str(n) + ["s", "p", "d", "f", "g", "h", "i"][l] + str(f))
        for n, l, f in config_tuples
    ]
    return " ".join(config_str)


def check_valid_quantum_number(Z, n, ell):
    """
    check the quantum number is valid
    Returns:
        valid or error
    """
    symbol = chemical_symbols[Z]
    config_tuple = config_str_to_config_tuples(load_orbital_configuration()[symbol])

    if not any([shell[:2] == (n, ell) for shell in config_tuple]):
        raise RuntimeError(
            f"Quantum numbers (n, ell) = ({n}, {ell}) not valid for element {symbol}"
        )