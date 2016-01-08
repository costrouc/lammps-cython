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
    cdef LAMMPS *_lammps
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

        self._lammps = new LAMMPS(argc, argv, comm)
        self.box = Box(self)
        self.atoms = Atoms(self)
        self.update = Update(self)
        self.thermo = Thermo(self)

    def __dealloc__(self):
        del self._lammps
        # TODO should I free string?

    @property
    def __version__(self):
        return self._lammps.universe.version

    def command(self, cmd):
        """Runs a LAMMPS command

           command: str
        """
        self._lammps.input.one(cmd)

    def file(self, filename):
        """Runs a LAMMPS file"""
        self._lammps.input.file(filename)

    def reset(self):
        self._lammps.input.one(b'clear')


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
    cdef UPDATE *_update
    def __cinit__(self, Lammps lammps):
        self._update = lammps._lammps.update

    @property
    def dt(self):
        return self._update.dt

    @property
    def time_step(self):
        return self._update.ntimestep

    @property
    def time(self):
        return self._update.atime


cdef class Atoms:
    cdef ATOM* _atom
    def __cinit__(self, Lammps lammps):
        self._atom = lammps._lammps._atom

    @property
    def num_total(self):
        return self._atom.natoms

    @property
    def num_local(self):
        return self._atom.nlocal

    def __len__(self):
        return self.num_local

    @property
    def tags(self):
        if self._atom.x == NULL:
            return None
        
        cdef size_t N = self.num_local
        cdef tagint[::1] array = <tagint[:N]>self._atom.tag
        return np.asarray(array)

    @property
    def positions(self):
        if self._atom.x == NULL:
            return None
        
        cdef size_t N = self.num_local
        cdef double[:, ::1] array = <double[:N, :3]>self._atom.x[0]
        return np.asarray(array)

    @property
    def velocities(self):
        if self._atom.v == NULL:
            return None

        cdef size_t N = self.num_local
        cdef double[:, ::1] array = <double[:N, :3]>self._atom.v[0]
        return np.asarray(array)

    @property
    def forces(self):
        if self._atom.f == NULL:
            return None
        
        cdef size_t N = self.num_local
        cdef double[:, ::1] arr = <double[:N, :3]>self._atom.f[0]
        return np.asarray(arr)

    @property
    def charges(self):
        if self._atom.q == NULL:
            return None
        
        cdef size_t N = self.num_local
        cdef double[::1] vector = <double[:N]>self._atom.q
        return np.asarray(vector)


cdef class Box:
    cdef DOMAIN* _domain
    def __cinit__(self, Lammps lammps):
        self._domain = lammps._lammps.domain

    @property
    def dimension(self):
        """ The dimension of the lammps run """
        return self._domain.dimension

    @property
    def lohi(self):
        cdef int dim = self.dimension
        cdef double[::1] boxlo = <double[:dim]>self._domain.boxlo
        cdef double[::1] boxhi = <double[:dim]>self._domain.boxhi
        return {'boxlo': np.array(boxlo), 'boxhi': np.array(boxhi)} # We copy arrays

    @property
    def tilts(self): #TODO how to handle 2d?
        cdef double xy = self._domain.xy
        cdef double xz = self._domain.xz
        cdef double yz = self._domain.yz
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
        cdef double vol = self._domain.xprd * self._domain.yprd
        if self.dimension == 2:
            return vol
        else: # dimension == 3
            return vol * self._domain.zprd
