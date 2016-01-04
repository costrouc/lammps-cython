#!/usr/bin/python3.4

import sys
sys.path.append('../../')

from lammps import Lammps

lmp = Lammps()

print("Version of Lammps: {}".format(lmp.__version__))
print("Number of atoms: {}".format(lmp.atoms.total_num))
print("Local number of atoms: {}".format(len(lmp.atoms)))
boxlo, boxhi = lmp.box.lengths
print(("Box Dimension:\n"
       "x {:3.2f} - {:3.2f}\n"
       "y {:3.2f} - {:3.2f}\n"
       "z {:3.2f} - {:3.2f}\n"
).format(boxlo[0], boxhi[0], boxlo[1], boxhi[1], boxlo[2], boxhi[2]))

lmp.file(b'in.melt')

print("Number of atoms: {}".format(lmp.atoms.total_num))
print("Local number of atoms: {}".format(len(lmp.atoms)))
boxlo, boxhi = lmp.box.lengths
print(("Box Dimension:\n"
       "x {:3.2f} - {:3.2f}\n"
       "y {:3.2f} - {:3.2f}\n"
       "z {:3.2f} - {:3.2f}\n"
).format(boxlo[0], boxhi[0], boxlo[1], boxhi[1], boxlo[2], boxhi[2]))

print(len(lmp.atoms.tags))
