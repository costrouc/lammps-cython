#!/usr/bin/python3.4

from lammps import Lammps

import matplotlib.pyplot as plt
import numpy as np

# lammps command list option'-sc', 'none' used so 
# lammps doesn't print tons of information
# lmp = Lammps(args=['lmp', '-sc', 'none'])
lmp = Lammps()

# print important initial lammps information
print((
    "Version of Lammps: {}\n"
    "Total number of atoms in simulation: {}\n"
    "Local number of atoms: {}\n"
).format(lmp.__version__, lmp.system.total, lmp.system.local))

# Read in lammps input file
lmp.file('in.melt')

print((
    "Total number of atoms: {}\n"
    "Local number of atoms: {}\n"
    "Time step:             {}\n"
).format(lmp.system.total, lmp.system.local, lmp.dt))

# Run LAMMPS simluation for 1000 steps
lmp.run(1000)

vel = lmp.system.velocities
plt.hist(np.linalg.norm(vel, axis=1), bins=40)
plt.savefig("vel.png")
