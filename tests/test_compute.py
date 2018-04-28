import functools

import lammps


# Redefine Lammps command-line args so no annoying logs or stdout
Lammps = functools.partial(lammps.Lammps, args=[
    '-log', 'none',
    '-screen', 'none'
])


def test_compute():
    lammps = Lammps()
