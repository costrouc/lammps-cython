#!/bin/bash

cat > lammps.cfg <<EOF
[lammps]
lammps_include_dir = $PREFIX/include/lammps
lammps_library_dir = $PREFIX/lib
lammps_library = lammps

[mpi]
mpi_include_dir = $PREFIX/include
mpi_library_dir = $PREFIX/lib
mpi_library = mpi, mpi_cxx

EOF

$PYTHON setup.py install
