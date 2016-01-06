include "lammps.pyd"

from libc.stdlib cimport malloc, free
cimport numpy as np

# Import the Python-level symbols
from mpi4py import MPI
import numpy as np

# Ownership of numpy array (2 steps)
# cdef extern from "numpy/arrayobject.h":
#     void PyArray_ENABLEFLAGS(np.ndarray arr, int flags)
#
# PyArray_ENABLEFLAGS(a, np.NPY_OWNDATA)

cdef class Lammps:
    cdef LAMMPS *thisptr
    cdef public Box box
    cdef public Atoms atoms
    def __cinit__(self):
        args = [] # TODO reimplement properly
        cdef int argc = len(args)
        cdef char** argv = <char**>malloc(len(args) * sizeof(char*))
        cdef mpi.MPI_Comm comm = mpi.MPI_COMM_WORLD
        for i in range(len(args)):
            temp_string = args[i].encode('UTF-8')
            argv[i] = temp_string # To prevent garbage collection

        self.thisptr = new LAMMPS(argc, argv, comm)
        self.box = Box(self)
        self.atoms = Atoms(self)

    def __dealloc__(self):
        del self.thisptr
        # TODO should I free string?

    @property
    def __version__(self):
        return self.thisptr.universe.version

    def command(self, cmd):
        """Runs a LAMMPS command

           command: str
        """
        self.thisptr.input.one(cmd)

    def file(self, filename):
        """Runs a LAMMPS file"""
        self.thisptr.input.file(filename)

    def reset(self):
        self.thisptr.input.one(b'clear')


cdef class Atoms:
    cdef LAMMPS *thisptr
    def __cinit__(self, Lammps lammps):
        self.thisptr = lammps.thisptr

    @property
    def total_num(self):
        return self.thisptr.atom.natoms

    @property
    def local_num(self):
        return self.thisptr.atom.nlocal

    def __len__(self):
        return self.local_num

    @property
    def tags(self):
        if self.thisptr.atom.x == NULL:
            return None
        
        cdef size_t N = self.local_num
        cdef tagint[::1] array = <tagint[:N]>self.thisptr.atom.tag
        return np.asarray(array)

    # TODO how to cast a double** array (why did lammps use this data structure)
    # Unique way that lammps represents internal data structure
    @property
    def positions(self):
        if self.thisptr.atom.x == NULL:
            return None
        
        cdef size_t N = self.local_num
        cdef double[:, ::1] array = <double[:N, :3]>self.thisptr.atom.x[0]
        return np.asarray(array)

    @property
    def velocities(self):
        if self.thisptr.atom.v == NULL:
            return None
        
        cdef size_t N = self.local_num
        cdef double[:, ::1] array = <double[:N, :3]>self.thisptr.atom.v[0]
        return np.asarray(array)

    @property
    def forces(self):
        if self.thisptr.atom.f == NULL:
            return None
        
        cdef size_t N = self.local_num
        cdef double[:, ::1] array = <double[:N, :3]>self.thisptr.atom.f[0]
        return np.asarray(array)

    @property
    def charges(self):
        if self.thisptr.atom.q == NULL:
            return None
        
        cdef size_t N = self.local_num
        cdef double[::1] array = <double[:N]>self.thisptr.atom.q
        return np.asarray(array)


cdef class Box:
    cdef LAMMPS *thisptr
    def __cinit__(self, Lammps lammps):
        self.thisptr = lammps.thisptr

    @property
    def dimension(self):
        """ The dimension of the lammps run """
        return self.thisptr.domain.dimension

    @property
    def lengths(self):
        cdef double[:] boxlo = self.thisptr.domain.boxlo
        cdef double[:] boxhi = self.thisptr.domain.boxhi
        return boxlo, boxhi

    @property
    def tilts(self):
        cdef double xy = self.thisptr.domain.xy
        cdef double xz = self.thisptr.domain.xz
        cdef double yz = self.thisptr.domain.yz
        return xy, xz, yz
