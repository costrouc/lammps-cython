from libc.stdlib cimport malloc, free
from libc.stdint cimport int64_t

from mpi4py import MPI          # Automatically calls Init and Finalize
cimport mpi4py.libmpi as mpi

# bigint is of type int64_t
# TODO do we need the constructors (excluding LAMMPS)?


cdef extern from "input.h" namespace "LAMMPS_NS":
     cdef cppclass Input:
          Input(LAMMPS*, int, char**) except +
          char* one(const char*)
          void file(const char*)

# cdef extern from "atom.h" namespace "LAMMPS_NS":
#      cdef cppclass Atom:
#           Atom(LAMMPS*) except +
#           int64_t natoms

cdef extern from "universe.h" namespace "LAMMPS_NS":
     cdef cppclass Universe:
          Universe(LAMMPS*, mpi.MPI_COMM) except +
          const char* version
          
# cdef extern from "domain.h" namespace "LAMMPS_NS":
#      cdef cppclass Domain:
#           Domain(LAMMPS*) except +
#           int box_exist
#           int dimension
#           int nonperiodic

cdef extern from "lammps.h" namespace "LAMMPS_NS":
    cdef cppclass LAMMPS:
         LAMMPS(int, char**, mpi.MPI_Comm) except +
         Input* input
         Universe* universe
         

cdef class Lammps:
    cdef LAMMPS *thisptr
    def __cinit__(self):
        args = [] # TODO reimplement properly
        cdef int argc = len(args)
        cdef char** argv = <char**>malloc(len(args) * sizeof(char*))
        cdef mpi.MPI_Comm comm = mpi.MPI_COMM_WORLD
        for i in range(len(args)):
            temp_string = args[i].encode('UTF-8')
            argv[i] = temp_string # To prevent garbage collection

        self.thisptr = new LAMMPS(argc, argv, comm)

    @property
    def __version__(self):
        return self.thisptr.universe.version

    def command(self, cmd):
        """Runs a LAMMPS command

           command: str
        """
        self.thisptr.input.one(cmd)

    def __dealloc__(self):
        del self.thisptr
        # TODO should I free string?
        
