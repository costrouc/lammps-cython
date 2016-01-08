#!/usr/bin/python3.4

import sys
sys.path.append('../../')

from lammps import Lammps

lmp = Lammps(['lmp', '-sc', 'none'])

print("Version of Lammps: {}".format(lmp.__version__))
print("Number of atoms: {}".format(lmp.atoms.num_total))
print("Local number of atoms: {}".format(len(lmp.atoms)))

lmp.file(b'in.melt')

print("Number of atoms: {}".format(lmp.atoms.num_total))
print("Local number of atoms: {}".format(len(lmp.atoms)))
