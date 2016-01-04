from libc.stdlib cimport malloc, free

from mpi4py import MPI          # Automatically calls Init and Finalize
cimport mpi4py.libmpi as mpi

cdef extern from "input.h" namespace "LAMMPS_NS":
     cdef cppclass Input:
          Input(LAMMPS*, int, char**) except +
          char* one(const char*)
          void file(const char*)

cdef extern from "lammps.h" namespace "LAMMPS_NS":
    cdef cppclass LAMMPS:
         LAMMPS(int, char**, mpi.MPI_Comm) except +
         Input* input
         int num_package
         

cdef class Lammps:
    cdef LAMMPS *thisptr
    def __cinit__(self, args):
        cdef int argc = len(args)
        cdef char** argv = <char**>malloc(len(args) * sizeof(char*))
        cdef mpi.MPI_Comm comm = mpi.MPI_COMM_WORLD
        for i in range(len(args)):
            temp_string = args[i].encode('UTF-8')
            argv[i] = temp_string # To prevent garbage collection

        self.thisptr = new LAMMPS(argc, argv, comm)

    def test(self):
        self.thisptr.input.one('print "Test that LAMMPS works!"')

    def __dealloc__(self):
        del self.thisptr
        # TODO should I free string?

    def get_packages(self):
        return self.thisptr.num_package
        
