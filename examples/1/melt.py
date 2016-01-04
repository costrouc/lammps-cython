import sys
sys.path.append('../../')

from lammps import Lammps

lmp = Lammps()

print("Version of Lammps: {}".format(lmp.__version__))
boxlo, boxhi = lmp.box.lengths
print(("Box Dimension:\n"
       "x {:3.2f} - {:3.2f}\n"
       "y {:3.2f} - {:3.2f}\n"
       "z {:3.2f} - {:3.2f}\n"
).format(boxlo[0], boxhi[0], boxlo[1], boxhi[1], boxlo[2], boxhi[2]))

lmp.file(b'in.melt')

boxlo, boxhi = lmp.box.lengths
print(("Box Dimension:\n"
       "x {:3.2f} - {:3.2f}\n"
       "y {:3.2f} - {:3.2f}\n"
       "z {:3.2f} - {:3.2f}\n"
).format(boxlo[0], boxhi[0], boxlo[1], boxhi[1], boxlo[2], boxhi[2]))
