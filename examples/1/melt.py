#!/usr/bin/python3.4

from lammps import Lammps

import matplotlib.pyplot as plt
import numpy as np

# lammps command list option'-sc', 'none' used so 
# lammps doesn't print tons of information
lmp = Lammps(args=['lmp', '-sc', 'none'])
# lmp = Lammps()

# print important initial lammps information
print(("Version of Lammps: {}\n"
       "Total number of atoms in simulation: {}\n"
       "Local number of atoms: {}\n").format(
           lmp.__version__,
           lmp.system.total,
           lmp.system.local))

# Read in lammps input file
lmp.file(b'in.melt')

print(("Total number of atoms: {}\n"
       "Local number of atoms: {}\n").format(
           lmp.system.total,
           lmp.system.local))

print(lmp.system.tags)
print(lmp.system.positions.shape)

for atom in lmp.system:
    print(atom.position)

# temps = []
# for i in range(100):
#     print(i)
#     lmp.run(1)
#     if lmp.thermo.temperature.scalar < 7.0:
#         print("ERROR: temp too low!!!")
#         exit()

#     temps.append(lmp.thermo.temperature.scalar)

# json.dump(temps, open("temp.dump", "w"))

# plt.plot(temps)
# plt.savefig("temp.png")

# # Do some interesting plotting
# vel = lmp.system.velocities

# plt.hist(np.linalg.norm(vel, axis=1), bins=40)
# plt.savefig("vel.png")
