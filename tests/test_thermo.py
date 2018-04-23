import os
import functools

import pytest
import lammps

# Redefine Lammps command-line args so no annoying logs or stdout
Lammps = functools.partial(lammps.Lammps, args=[
    '-log', 'none',
    '-screen', 'none'
])


def test_thermo_init_computes():
    lmp = Lammps()

    assert len(lmp.thermo.computes) == 3

    assert lmp.thermo.computes['thermo_temp'].style == 'temp'
    assert lmp.thermo.temperature.style == 'temp'
    assert lmp.thermo.temperature.name == 'thermo_temp'

    assert lmp.thermo.computes['thermo_press'].style == 'pressure'
    assert lmp.thermo.pressure.style == 'pressure'
    assert lmp.thermo.pressure.name == 'thermo_press'

    assert lmp.thermo.computes['thermo_pe'].style == 'pe'
    assert lmp.thermo.potential_energy.style == 'pe'
    assert lmp.thermo.potential_energy.name == 'thermo_pe'


def test_thermo_add_compute():
    lmp = Lammps()

    # com (center of mass)
    thermo_id = 'thermo_test'
    thermo_style = 'temp'
    lmp.thermo.add(id=thermo_id, style=thermo_style)

    assert lmp.thermo.computes[thermo_id].name == thermo_id
    assert lmp.thermo.computes[thermo_id].style == thermo_style
