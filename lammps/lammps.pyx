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
    cdef public Update update
    cdef public Thermo thermo
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
        self.update = Update(self)
        self.thermo = Thermo(self)

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


cdef class Thermo:
    cdef LAMMPS *thisptr
    cdef public Compute temperature
    cdef public Compute pressure
    cdef public Compute energy
    def __cinit__(self, Lammps lammps):
        self.thisptr = lammps.thisptr

        cdef int index_temp = self.thisptr.modify.find_compute(b"thermo_temp")
        cdef int index_press = self.thisptr.modify.find_compute(b"thermo_press")
        cdef int index_pe = self.thisptr.modify.find_compute(b"thermo_pe")

        self.temperature = Compute(lammps, index_temp)
        self.pressure = Compute(lammps, index_press)
        self.energy = Compute(lammps, index_pe)



cdef class Compute:
     cdef COMPUTE *thisptr
     def __cinit__(self, Lammps lammps, int index):
         self.thisptr = lammps.thisptr.modify.compute[index]

     @property
     def scalar(self):
         if self.thisptr.scalar_flag == 0:
             raise NotImplementedError()

         self.thisptr.compute_scalar()
         return self.thisptr.scalar

     @property
     def vector(self):
         if self.thisptr.vector_flag == 0:
             raise NotImplementedError()

         cdef int n = self.thisptr.size_vector
         self.thisptr.compute_vector()
         cdef double[::1] array = <double[:n]>self.thisptr.vector
         return np.asarray(array)
         

cdef class Update:
    cdef LAMMPS *thisptr
    def __cinit__(self, Lammps lammps):
        self.thisptr = lammps.thisptr

    @property
    def dt(self):
        return self.thisptr.update.dt

    @property
    def time_step(self):
        return self.thisptr.update.ntimestep

    @property #TODO not exact (see update.atimestep)
    def time(self):
        return self.thisptr.update.atime


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
    def lohi(self):
        cdef int dim = self.dimension
        cdef double[::1] boxlo = <double[:dim]>self.thisptr.domain.boxlo
        cdef double[::1] boxhi = <double[:dim]>self.thisptr.domain.boxhi
        return {'boxlo': np.array(boxlo), 'boxhi': np.array(boxhi)} # We copy arrays

    @property
    def tilts(self): #TODO how to handle 2d?
        cdef double xy = self.thisptr.domain.xy
        cdef double xz = self.thisptr.domain.xz
        cdef double yz = self.thisptr.domain.yz
        return {'xy': xy, 'xz': xz, 'yz': yz}
    
    # See http://lammps.sandia.gov/doc/Section_howto.html#4_12
    @property
    def lengths(self):
        raise NotImplementedError()

    @property
    def angles(self):
        raise NotImplementedError()

    @property
    def volume(self):
        cdef double vol = self.thisptr.domain.xprd * self.thisptr.domain.yprd
        if self.dimension == 2:
            return vol
        else: # dimension == 3
            return vol * self.thisptr.domain.zprd
