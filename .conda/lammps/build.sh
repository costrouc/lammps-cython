#!/bin/bash

export PACKAGES="asphere body class2 colloid compress coreshell dipole granular kspace manybody molecule mc misc opt peri qeq replica rigid shock snap srd user-reaxc"
export BUILD="mpi"
export LMP_INCLUDES="-DLAMMPS_EXCEPTIONS -DLAMMPS_GZIP -DLAMMPS_MEMALIGN=64"
export NUM_CPUS=2

# Make binary, shared library, and static library
cd src
for pack in $PACKAGES; do make "yes-$pack"; done
make mode=lib $BUILD -j${NUM_CPUS} LMP_INC="$LMP_INCLUDES"
make mode=shlib $BUILD -j${NUM_CPUS} LMP_INC="$LMP_INCLUDES"
make $BUILD -j${NUM_CPUS} LMP_INC="$LMP_INCLUDES"

# copy binary, shared library, static library, and include files
cp lmp_$BUILD $PREFIX/bin/lammps
cp liblammps_$BUILD.so $PREFIX/lib/liblammps.so
cp liblammps_$BUILD.a $PREFIX/lib/liblammps.a
mkdir -p $PREFIX/include/lammps
cp *.h $PREFIX/include/lammps/
