"""LAMMPS Python interface (the GOOD parts) 

This interface is inspired by HOOMD and tries its best to look
similar. I try to include only orthogonal functions (e.g. make
there only be ONE way to do something).

.. todo:: Features
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
    
..  todo:: determine if I need to free char* strings
    """
    cdef char** argv = <char**>malloc(len(args) * sizeof(char*))
    cdef int i
    for i in range(len(args)):
        temp = args[i].encode('UTF-8')
        argv[i] = temp
    return argv


cdef class Lammps:
    """LAMMPS base class

..  py:function:: __init__(args, comm=None)
    
    Initialize a Lammps object. 

    :param list args: list of command line args that would be supplied to normal lammps executable
    :param comm: mpi4py comm object default value is MPI_COMM_WORLD

    To see a list of possible lammps args (e.g. `command line
    arguments
    <http://lammps.sandia.gov/doc/Section_start.html#command-line-options>`_). These
    must be provided as a list of strings.
    """
    cdef LAMMPS *_lammps
    cdef mpi.MPI_Comm _comm
    cdef public Box box
    cdef public System system
    cdef public Thermo thermo
    def __cinit__(self, args, comm=None):
        """ Docstring in Lammps base class (sphinx can find doc when compiled) """
        cdef int argc = len(args)
        cdef char** argv = <char**>args_to_cargv(args)

        if comm is None:
            self._comm = mpi.MPI_COMM_WORLD
        else:
            raise NotImplementedError()

        self._lammps = new LAMMPS(argc, argv, self._comm)
        # don't free char** becuase used by lammps (they don't copy their strings!)

        self.box = Box(self)
        self.system = System(self)
        self.thermo = Thermo(self)

    def __dealloc__(self):
        del self._lammps

    @property
    def __version__(self):
        """ Prints the version of LAMMPS 

        Format is <day><month><year> e.g. 7Dec15
        """
        return self._lammps.universe.version

    def command(self, cmd):
        """Runs any single LAMMPS command

           :param str command: command for lammps to execute

        See lammps documentation for available `commands <http://lammps.sandia.gov/doc/Section_commands.html#individual-commands>`_.
        """
        self._lammps.input.one(cmd)

    def file(self, filename):
        """ Runs a LAMMPS input file

        :param str filename: filename of input file

        This is equivalent to setting the -i command line flag.

..      todo:: learn how file() behaves with multiple invocations
        """
        self._lammps.input.file(filename)

    def run(self, long steps):
        """ Runs the lammps simulation for N steps

        :param int steps: number of steps to run simulation equivalent to "run <steps>" command

        See lammps documentation for description of `run command
        <http://lammps.sandia.gov/doc/run.html>`_
        """
        self._lammps.input.one(str.encode('run {}'.format(steps)))

    def reset(self):
        """ Resets the lammps simulation

        Deletes all atoms, restores all settings to their default
        values, and frees all memory allocated by LAMMPS. Equivalent
        to the LAMMPS `clear command
        <http://lammps.sandia.gov/doc/clear.html>`_.
        """
        self._lammps.input.one(b'clear')

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


cdef class Thermo:
    """ Compute thermodynamic properties of a group of particles in system from available computes. 

    You must first define a compute before you can extract
    thermodynamics properties of the current time step. Three computes
    are always created, named “thermo_temp”, “thermo_press”, and
    “thermo_pe” these are initialized in the output.cpp in LAMMPS.

..  py:function:: __init__(Lammps)
    
    Initialize a Thermo object.

    :param lammps: Lammps object

..  py:attribute:: computes
    
    A dictionary of computes. {id: :py:class:Compute}.
    """
    cdef MODIFY* _modify
    cdef public Compute temperature
    cdef public Compute pressure
    cdef public Compute potential_energy
    cdef public dict computes
    def __cinit__(self, Lammps lammps):
        """ Docstring in Thermo base class (sphinx can find doc when compiled) """
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

    def add(self, id, style, args):
        """ Add a compute to LAMMPS

        :param str id: name of new lammps compute cannot conflict with existing compute ids
        :param str style: name of compute to add 

        See `compute <http://lammps.sandia.gov/doc/compute.html>`_ for
        more information on creating computes.
        """
        raise NotImplementedError()


cdef class Compute:
    """ Extract compute object property from LAMMPS system

    See `compute <http://lammps.sandia.gov/doc/compute.html>`_ for
    more information on available computes.

..  py:function:: __init__(Lammps)
    
    Initialize a Compute object.

    :param lammps: Lammps object
    """
    cdef COMPUTE *_compute
    def __cinit__(self, Thermo thermo, id, args=None):
        """ Docstring in Compute base class (sphinx can find doc when compiled) """
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
        """ Modify LAMMPS compute
        
        :param list args: list[str] of args

        See `compute_modify
        <http://lammps.sandia.gov/doc/compute_modify.html>`_ for more
        information on args etc.
        """
        cdef char** argv = args_to_cargv(args)
        cdef int argc = len(args)
        self._compute.modify_params(argc, argv)

    @property
    def id(self):
        """ name of compute id 

        :return: compute id
        :rtype: str  
        """
        return self._compute.id

    @property
    def style(self):
        """ style of compute id 

        :return: compute style
        :rtype: str
        """
        return self._compute.style

    @property
    def scalar(self):
        """ scalar value of compute

        :return: value of compute
        :rtype: float
        :raises NotImplementedError: if LAMMPS does not have a scalar function for this compute
        """
        if self._compute.scalar_flag == 0:
            raise NotImplementedError("style {} does not have a scalar function".format(self.style))

        self._compute.compute_scalar()
        return self._compute.scalar

    @property
    def vector(self):
        """ vector value of compute

        :return: vector value of compute
        :rtype: numpy.ndarray
        :raises NotImplementedError: if LAMMPS does not have a vector function for this compute
        """
        if self._compute.vector_flag == 0:
            raise NotImplementedError("style {} does not have a vector function".format(self.style))

        self._compute.compute_vector()

        cdef int N = self._compute.size_vector
        cdef double[::1] vector = <double[:N]>self._compute.vector
        return np.asarray(vector)
         

cdef class System:
    cdef ATOM* _atom
    def __cinit__(self, Lammps lammps):
        """ Docstring in Atom base class (sphinx can find doc when compiled) """
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
        """ Docstring in Box base class (sphinx can find doc when compiled) """
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
