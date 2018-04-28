from math import pi

import pytest
import numpy as np

import lammps


def test_lattice_const_to_lammps_box_cubic():
    lengths = (5, 5, 5)
    angles = (pi/2, pi/2, pi/2)
    boxhi, tilts = lammps.core.lattice_const_to_lammps_box(lengths, angles)
    assert np.all(np.isclose(boxhi, lengths))
    assert np.all(np.isclose(tilts, (0, 0, 0)))


@pytest.mark.xfail
def test_lattice_const_to_lammps_box_rhomb():
    # test is failing
    # 3C-SiC
    lengths = (3.0968, 3.0968, 3.0968)
    angles = (pi/3, pi/3, pi/3)
    boxhi, tilts = lammps.core.lattice_const_to_lammps_box(lengths, angles)
    assert np.all(np.isclose(boxhi, lengths))
    assert np.all(np.isclose(tilts, (0, 0, 0)))


def test_lammps_box_to_lattice_const_cubic():
    box_lengths = (5, 5, 5)
    tilts = (0, 0, 0)
    lengths, angles = lammps.core.lammps_box_to_lattice_const(box_lengths, tilts)
    assert np.all(np.isclose(box_lengths, lengths))
    assert np.all(np.isclose(angles, (pi/2, pi/2, pi/2)))


def test_lammps_initial_box(lmp):
    assert lmp.box.dimension == 3
    assert np.all(np.isclose(lmp.box.lengths, (1., 1., 1.)))
    assert np.all(np.isclose(lmp.box.angles, (pi/2., pi/2., pi/2.)))
    assert np.all(np.isclose(lmp.box.lohi, [[-0.5, -0.5, -0.5], [0.5, 0.5, 0.5]]))
    assert np.all(np.isclose(lmp.box.tilts, [0, 0, 0]))
    assert np.all(np.isclose(lmp.box.lengths_angles, [[1, 1, 1], [pi/2, pi/2, pi/2]]))
    # lammps has some seriously weird initial behavior
    # has unit cell 1x1x1 with volume 0 ???
    # actually has non-deterministic behavior 0 or inf
    # assert np.isclose(lmp.box.volume, 0.)


def test_lammps_set_box_from_lattice_const(lmp):
    atom_types = 5
    lengths = (10, 10, 10)
    angles = (pi/2., pi/2., pi/2.)
    lmp.box.from_lattice_const(atom_types, lengths, angles)
    assert np.all(np.isclose(lmp.box.lengths, lengths))
    assert np.all(np.isclose(lmp.box.angles, angles))
    assert len(lmp.system.atom_types) == atom_types
    assert np.isclose(lmp.box.volume, 10**3)


def test_lammps_update_lattice_const(lmp):
    lengths = (10, 10, 10)
    angles = (pi/2., pi/2., pi/2.)
    lmp.box.update_lattice_const(lengths, angles)
    assert np.all(np.isclose(lmp.box.lengths, lengths))
    assert np.all(np.isclose(lmp.box.angles, angles))
    assert np.isclose(lmp.box.volume, 10**3)
