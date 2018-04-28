""" Test integration for groups needs

This allows me to ensure that particular workflows are working
"""

from math import pi

import numpy as np


def test_dr_xu_integration(lmp):
    """ Dr. Xu integration requirements

     - set unit cell
    - set num atoms
    - set atom mass
    - create atoms
    - run similation
    - global gather atoms properties (x, f, ...)
    - global set atom properties (x, f, ...)
    """
    lengths = (4.2, 4.2, 4.2)
    angles = (pi/2, pi/2, pi/2)
    elements = [{'symbol': 'Mg', 'mass': 24.3}, {'symbol': 'O', 'mass': 16.0}]
    symbol_indicies = {element['symbol']: i+1 for i, element in enumerate(elements)}
    symbols = ['Mg', 'Mg', 'Mg', 'Mg', 'O', 'O', 'O', 'O']

    positions = [
        (0, 0, 0), (2.1, 2.1, 0), (2.1, 0, 2.1), (0, 2.1, 2.1), # Mg
        (2.1, 0, 0), (0, 2.1, 0), (0, 0, 2.1), (2.1, 2.1, 2.1)  # O
    ]

    # Setup Unit cell
    lmp.box.from_lattice_const(len(elements), lengths, angles)
    assert np.all(np.isclose(lmp.box.lengths, lengths))
    assert np.all(np.isclose(lmp.box.angles, angles))
    assert len(lmp.system.atom_types) == len(elements)

    # Set Atom Masses
    assert np.all(np.isclose(np.array([a.mass for a in lmp.system.atom_types]), [0.] * len(elements)))
    for element, atom_type in zip(elements, lmp.system.atom_types):
        atom_type.mass = element['mass']
    assert np.all(np.isclose(np.array([a.mass for a in lmp.system.atom_types]), [e['mass'] for e in elements]))

    # Set atom positions and symbols
    assert lmp.system.total == 0
    atom_types = np.array([symbol_indicies[symbol] for symbol in symbols], dtype=np.int32)
    velocities = np.zeros((len(atom_types), 3))
    lmp.system.create_atoms(atom_types, np.array(positions), velocities)
    assert lmp.system.total == len(atom_types)
    # more tests on values
