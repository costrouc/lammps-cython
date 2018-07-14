import numpy as np


def write_table_pair_potential(func, dfunc=None, bounds=(1.0, 10.0), samples=1000, toll=1e-6, keyword='PAIR', filename='lammps.table'):
    """A helper function to write lammps pair potentials. Assumes that functions are vectorized.

    """
    r_min, r_max = bounds
    if dfunc is None:
        dfunc = lambda r: (func(r+toll) - func(r-toll)) / (2*toll)

    with open(filename, 'w') as f:
        i = np.arange(1, n+1)
        r = np.linspace(r_min, r_max, n)
        forces = func(r)
        energies = dfunc(r)
        lines = ['%d %f %f %f\n' % (index, radius, force, energy) for index, radius, force, energy in zip(i, r, forces, energies)]
        f.write("%s\nN %d\n\n" % (keyword, samples) + ''.join(lines))
