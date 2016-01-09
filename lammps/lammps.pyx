"""LAMMPS Python interface (the GOOD parts) 

This interface is inspired by HOOMD and tries its best to look
similar. I try to include only orthogonal functions (e.g. make
there only be ONE way to do something).

todo Features
 - Group support
 - integrator support
 - finish not implemented functions

"""


include "lammps.pyd"

from libc.stdlib cimport malloc, free
cimport numpy as np

# Import the Python-level symbols
from mpi4py import MPI
import numpy as np


# helper functions char**
cdef char** args_to_cargv(args):
    """ Convert list of args[str] to char** 
    
    TODO: determine if I need to free char* strings
    """
    cdef char** argv = <char**>malloc(len(args) * sizeof(char*))
    cdef int i
    for i in range(len(args)):
        temp = args[i].encode('UTF-8')
        argv[i] = temp
    return argv


cdef class Lammps:
    r""" LAMMPS base class

    Representation of LAMMPS
       box 

    :math:`a^2 + b^2 = c^2`
    :math:`\alpha \beta \frac{1}{2}`
    :math:`\alpha + \beta  = \frac{1}{2}`
    
..  todo:: Something to do.
    """
    cdef LAMMPS *_lammps
    cdef mpi.MPI_Comm _comm
    cdef public Box box
    cdef public Atoms system
    cdef public Update update
    cdef public Thermo thermo
    def __cinit__(self, args, comm=None):
        cdef int argc = len(args)
        cdef char** argv = <char**>args_to_cargv(args)

        if comm is None:
            self._comm = mpi.MPI_COMM_WORLD
        else:
            raise NotImplementedError()

        self._lammps = new LAMMPS(argc, argv, self._comm)
        # don't free char** becuase used by lammps (they don't copy their strings!)

        self.box = Box(self)
        self.system = Atoms(self)
        self.update = Update(self)
        self.thermo = Thermo(self)

    def __dealloc__(self):
        del self._lammps

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

    def run(self, long steps):
        self._lammps.input.one(str.encode('run {}'.format(steps)))

    def send_message(sender, recipient, message_body, priority):
        """ Send message (not true)
        
        :param str sender: The person sending the message
        :param str recipient: The recipient of the message
        :param str message_body: The body of the message
        :param priority: The priority of the message, can be a number 1-5
        :type priority: integer or None
        :return: the message id
        :rtype: int
        :raises ValueError: if the message_body exceeds 160 characters
        :raises TypeError: if the message_body is not a basestring
        """
        return 10

    def reset(self):
        self._lammps.input.one(b'clear')


cdef class Thermo:
    cdef MODIFY* _modify
    cdef public Compute temperature
    cdef public Compute pressure
    cdef public Compute potential_energy
    cdef public dict computes
    def __cinit__(self, Lammps lammps):
        self._modify = lammps._lammps.modify
        self.computes = dict()

        # Add computes automatically added by
        # Lammps (output.h)
        self.temperature = Compute(self, b"thermo_temp")
        self.computes.update({'thermo_temp': self.temperature})
        self.pressure = Compute(self, b"thermo_press")
        self.computes.update({'thermo_press': self.pressure})
        self.potential_energy = Compute(self, b"thermo_pe")
        self.computes.update({'thermo_pe': self.potential_energy})

    def modify(self, id, args):
        raise NotImplementedError()

    def add(self, id, args):
        raise NotImplementedError()


cdef class Compute:
     cdef COMPUTE *_compute
     def __cinit__(self, Thermo thermo, id, args=None):
         cdef int index = thermo._modify.find_compute(id)
         cdef char** argv
         cdef int argc
         if index == -1 and args is None:
             raise ValueError("args must be supplied for new compute")
         elif index != -1 and args:
             raise ValueError("compute id already exists use modify or delete")
         elif index == -1:
             # Hack to add id to args (needed for add compute)
             args = [id] + args
             argv = args_to_cargv(args)
             argc = len(args)
             thermo._modify.add_compute(argc, argv)
             index = thermo._modify.find_compute(id)
        
         self._compute = thermo._modify.compute[index]

     def modify(self, args):
         cdef char** argv = args_to_cargv(args)
         cdef int argc = len(args)
         self._compute.modify_params(argc, argv)

     @property
     def id(self):
         return self._compute.id

     @property
     def style(self):
         return self._compute.style

     @property
     def scalar(self):
         if self._compute.scalar_flag == 0:
             raise NotImplementedError()

         self._compute.compute_scalar()
         return self._compute.scalar

     @property
     def vector(self):
         if self._compute.vector_flag == 0:
             raise NotImplementedError()

         self._compute.compute_vector()

         cdef int N = self._compute.size_vector
         cdef double[::1] vector = <double[:N]>self._compute.vector
         return np.asarray(vector)
         

cdef class Update:
    cdef UPDATE *_update
    def __cinit__(self, Lammps lammps):
        self._update = lammps._lammps.update

    @property
    def dt(self):
        return self._update.dt

    @dt.setter
    def dt(self, double value):
        self._update.dt = value

    @property
    def time_step(self):
        return self._update.ntimestep

    @time_step.setter
    def time_step(self, bigint value):
        self._update.reset_timestep(value)
    
    @property
    def time(self):
        return self._update.atime


cdef class Atoms:
    cdef ATOM* _atom
    def __cinit__(self, Lammps lammps):
        self._atom = lammps._lammps.atom

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
