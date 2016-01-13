#!python
#cython: embedsignature=True
"""LAMMPS Python interface (the GOOD parts) 

This interface is inspired by HOOMD and tries its best to look
similar. I try to include only orthogonal functions (e.g. make
there only be ONE way to do something).

.. todo:: Features
 - Group support
 - integrator support
 - lammps units support (maybe integrate with units package)
 - lammps construct box from lengths and angles


Notes:
Creating from atoms and unit cell

atom_style (atomic or charge)

region lammps-python prism xlo xhi ylo yhi zlo zhi xy xz yz
create_box <number of atom types> <region-id>

# create_atoms <atom_type> single <x> <y> <z> 
#  - use remap to put atoms outside cell inside
#  - units [box or lattice] - fractional or real

No we use the atom_vec.create_atom(int, double*) call!!
called by atom->avec->create_atom
only create atom if it is in subbox of processor
"""

include "core.pyd"

from libc.stdlib cimport malloc, free
cimport numpy as np
from libc.math cimport sqrt, acos

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

..  py:function:: __init__(self, args, comm=None)
    
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
    def __cinit__(self, args=None, comm=None):
        """ Docstring in Lammps base class (sphinx can find doc when compiled) """
        if comm is None:
            comm = MPI.COMM_WORLD

        if not isinstance(comm, MPI.Comm):
            raise TypeError("comm arg must be of type MPI.Comm")

        # From discussion (todo: possible bug here)
        # https://groups.google.com/forum/#!topic/mpi4py/jPqNrr_8UWY
        # https://bitbucket.org/mpi4py/mpi4py/issues/13/import-failure-on-cython-extensions-that
        # https://bitbucket.org/mpi4py/mpi4py/commits/93804f5609ceac5e42aad58e998ebb1f213296c8
        cdef size_t comm_addr= MPI._addressof(comm)
        cdef mpi.MPI_Comm* comm_ptr = <mpi.MPI_Comm*>comm_addr
        self._comm = comm_ptr[0]
        
        if args is None:
            args = ['python']
        
        cdef int argc = len(args)
        cdef char** argv = <char**>args_to_cargv(args)

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
        """ timestep size for run step in simulation time units

        :getter: Returns the timestep size
        :setter: Sets the timestep size
        """
        return self._lammps.update.dt

    @dt.setter
    def dt(self, double value):
        self._lammps.update.dt = value

    @property
    def time_step(self):
        """ current number of timesteps that have been run
    
        :getter: Returns the timestep
        :setter: Sets the timestep
        """
        return self._lammps.update.ntimestep

    @time_step.setter
    def time_step(self, bigint value):
        self._lammps.update.reset_timestep(value)
    
    @property
    def time(self):
        """ total time that has elapsed from lammps runs in simulation time units
    
        :getter: Returns the total time
        :setter: Sets the total time
        """
        return self._lammps.update.atime


cdef class Thermo:
    """ Compute thermodynamic properties of a group of particles in system from available computes. 

    You must first define a compute before you can extract
    thermodynamics properties of the current time step. Three computes
    are always created, named “thermo_temp”, “thermo_press”, and
    “thermo_pe” these are initialized in the output.cpp in LAMMPS.

..  py:function:: __init__(self, Lammps)
    
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

..  py:function:: __init__(self, Lammps)
    
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


cdef class Atom:
    """ Represents a single atom in the LAMMPS simulation

    Since LAMMPS is a distributed system each atom only exists on 
    one of the processors.

..  py:function:: __init__(self, System)

    :param int tag: tag number of particle 
    :param int lindex: the local index of the atom

    Initialize an System object. Must supply either a tag number or index
    """
    cdef ATOM* _atom
    cdef readonly long tag
    cdef long local_index
    def __cinit__(self, System system, tagint tag=0, lindex=None):
        """ Docstring in System base class (sphinx can find doc when compiled) """
        self._atom = system._atom

        if lindex is not None:
            if lindex < 0 or lindex >= system.local:
                raise IndexError("index not withing local atoms index")
            self.local_index = lindex
            self.tag = self._atom.tag[self.local_index]
        else:
            raise NotImplementedError("currently need the local index")
            # self.tag = tag
            # self.index = self._atom.map(self.tag)

    @property
    def type(self):
        """ type of local atom
        
        :getter: returns the type of the local atom
        :setter: sets the type of the local atom
        """
        return self._atom.type[self.local_index]

    @type.setter
    def type(self, value):
        self._atom.type[self.local_index] = value

    @property
    def position(self):
        """ position of local atom
        
        :getter: returns the position of the local atom
        """
        # Check if atom is on this processor
        if (self.local_index == -1) or (self._atom.x == NULL):
            return None

        cdef double[::1] array = <double[:3]>&self._atom.x[0][self.local_index * 3]
        return np.asarray(array)

    @property
    def velocity(self):
        """ velocity of local atom
        
        :getter: returns the velocity of the local atom
        """
        # Check if atom is on this processor
        if (self.local_index == -1) or (self._atom.v == NULL):
            return None

        cdef double[::1] array = <double[:3]>&self._atom.v[0][self.local_index * 3]
        return np.asarray(array)

    @property
    def force(self):
        """ force of local atom
        
        :getter: returns the force of the local atom
        """
        # Check if atom is on this processor or property is set
        if (self.local_index == -1) or (self._atom.f == NULL):
            return None

        cdef double[::1] array = <double[:3]>&self._atom.f[0][self.local_index * 3]
        return np.asarray(array)

    @property
    def charge(self):
        """ charge of local atom
        
        :getter: returns the charge of the local atom
        :setter: sets the chage of the local atom
        """
        # Check if atom is on this processor or property is set
        if (self.local_index == -1) or (self._atom.q == NULL):
            return None

        return self._atom.q[self.local_index]

    @charge.setter
    def charge(self, value):
        # Check if atom is on this processor or property is set
        if (self.local_index == -1) or (self._atom.q == NULL):
            return None

        self._atom.q[self.local_index] = value
         

cdef class System:
    """ Represents all the atoms in the LAMMPS simulation

    Since LAMMPS is a distributed system each processor has a local
    view of its Atoms.

..  py:function:: __init__(self, Lammps)
    
    Initialize a System object.

    :param lammps: Lammps object
    """
    cdef ATOM* _atom
    cdef Lammps lammps
    cdef unsigned int local_index # used by iter
    def __cinit__(self, Lammps lammps):
        """ Docstring in System base class (sphinx can find doc when compiled) """
        self._atom = lammps._lammps.atom
        self.lammps = lammps
        self.local_index = 0

    @property
    def total(self):
        """ Total number of atoms in LAMMPS simulation

        :getter: Returns the total number of atoms in LAMMPS simulation
        """
        return self._atom.natoms

    @property
    def local(self):
        """ Local number of atoms stored in core

        :getter: Returns the local number of atoms specific to core
        """
        return self._atom.nlocal

    def __len__(self):
        return self.local

    def __getitem__(self, lindex):
        return Atom(self, lindex=lindex)

    def __iter__(self):
        self.local_index = 0
        return self

    def __next__(self):
        if self.local_index >= len(self):
            raise StopIteration()
        
        atom = Atom(self, self._atom.tag[self.local_index], lindex=self.local_index)
        self.local_index += 1
        return atom

    @property
    def tags(self):
        """ Tags associated with local atoms stored on core in numpy.ndarray

        :getter: Returns the local tags of atoms specific to core
        """
        if self._atom.x == NULL:
            return None
        
        cdef size_t N = self.local
        cdef tagint[::1] array = <tagint[:N]>self._atom.tag
        return np.asarray(array)

    @property
    def types(self):
        """ Tags associated with local atoms stored on core in numpy.ndarray

        :getter: Returns the local tags of atoms specific to core
        """
        if self._atom.x == NULL:
            return None
        
        cdef size_t N = self.local
        cdef int[::1] array = <int[:N]>self._atom.type
        return np.asarray(array)

    @property
    def positions(self):
        """ Positions associated with local atoms stored on core in numpy.ndarray

        :getter: Returns the local positions of atoms specific to core
        """
        if self._atom.x == NULL:
            return None
        
        cdef size_t N = self.local
        cdef double[:, ::1] array = <double[:N, :3]>self._atom.x[0]
        return np.asarray(array)

    @property
    def velocities(self):
        """ Velocities associated with local atoms stored on core in numpy.ndarray

        :getter: Returns the local velocities of atoms specific to core
        """
        if self._atom.v == NULL:
            return None

        cdef size_t N = self.local
        cdef double[:, ::1] array = <double[:N, :3]>self._atom.v[0]
        return np.asarray(array)

    @property
    def forces(self):
        """ Forces associated with local atoms stored on core

        :getter: Returns the local forces of atoms specific to core
        """
        if self._atom.f == NULL:
            return None
        
        cdef size_t N = self.local
        cdef double[:, ::1] arr = <double[:N, :3]>self._atom.f[0]
        return np.asarray(arr)

    @property
    def charges(self):
        """ Charges associated with local atoms stored on core

        :getter: Returns the local charges of atoms specific to core
        """
        if self._atom.q == NULL:
            return None

        cdef size_t N = self.local
        cdef double[::1] vector = <double[:N]>self._atom.q
        return np.asarray(vector)

    def add(self, atom_type, position, remap=False, units='box'):
        """Create atom in system

        :param int atom_type: atom type id
        :param np.array[3] position: position of atom
        :param boolean remap: whether to remap atoms in box or not
        :param str units: units to use for positions

        essentially runs command:
        create_atoms <atom_type> single <x> <y> <z>

        The lammps internals are hard to use. Yes there are some
        expensive mpi calls becuase of this (so this is slow for large
        systems). HACK

        See lammps `create_atoms
        <http://lammps.sandia.gov/doc/create_atoms>`_. documentation.
        """
        if len(position) != 3:
            raise ValueError('position array must be 3 values')

        cmd = (
            'create_atoms {!s} single {!s} {!s} {!s} units {!s} remap {!s}'
        ).format(atom_type, position[0], position[1], position[2], units,
                 'yes' if remap else 'no')
        self.lammps.command(cmd.encode('utf-8'))


cdef class Box:
    """ Represents the shape of the simulation cell.

..  py:function:: __init__(self, Lammps)
    
    Initialize a Box object.

    :param lammps: Lammps object
    """
    cdef DOMAIN* _domain
    def __cinit__(self, Lammps lammps):
        """ Docstring in Box base class (sphinx can find doc when compiled) """
        self._domain = lammps._lammps.domain

    @property
    def dimension(self):
        """ The dimension of the lammps run either 2D or 3D

        :returns: dimension of lammps simulation
        :rtype: int
        """
        return self._domain.dimension

    @property
    def lohi(self):
        r""" LAMMPS box description of boxhi and boxlo

        :return: dictionary of lower and upper in each dimension
        :rtype: dict

        For example a return dictionary would be

..      code-block:: python
        
        lohi = {
            'boxlo': np.array([0.0, 0.0, 0.0]),
            'boxhi': np.array([10.0, 10.0, 10.0])
        }

..      math::
        lx &= boxhi[0] - boxlo[0] = a \\
        ly &= boxhi[1] - boxlo[1] = \sqrt{b^2 - xy^2} \\
        lz &= boxhi[2] - boxlo[2] = \sqrt{c^2 - xz^2 - yz^2}


Additional documentation can be found at `lammps <http://lammps.sandia.gov/doc/Section_howto.html?highlight=prism#howto-12>`_.
        """
        cdef int dim = self.dimension
        cdef double[::1] boxlo = <double[:dim]>self._domain.boxlo
        cdef double[::1] boxhi = <double[:dim]>self._domain.boxhi
        return {'boxlo': np.array(boxlo), 'boxhi': np.array(boxhi)} # We copy arrays

    @property
    def tilts(self):
        r""" LAMMPS box description of xy xz yz

        :return: dictionary of xy xz yz
        :rtype: dict

        For example a return dictionary would be

..      code-block:: python

        tilts = {
           'xy': 0.0,
           'xz': 0.0,
           'yz': 0.0
        }

..      math::
        xy &= b \cos{\gamma} \\
        xz &= c \cos{\beta}  \\
        yz &= \frac{b c \cos{\alpha} - xy xz}{ly}


Additional documentation can be found at `lammps <http://lammps.sandia.gov/doc/Section_howto.html?highlight=prism#howto-12>`_.

..      todo::
         how to handle 2d? Needs to be addressed throughout
        """
        cdef double xy = self._domain.xy
        cdef double xz = self._domain.xz
        cdef double yz = self._domain.yz
        return {'xy': xy, 'xz': xz, 'yz': yz}

    @property
    def lengths_angles(self):
        r""" Calculates the lengths and angles from lammps box

        :getter: return lengths (a, b, c) and angles (alpha, beta, gamma)

        For some reason lammps decided to use a non standard box
        definition

..      math::

        a &= lx \\
        b &= \sqrt{ly^2 + xy^2} \\
        c &= \sqrt{lz^2 + xz^2 + yz^2} \\
        \cos{\alpha} &= \frac{xy xz + ly yz}{b c} \\
        \cos{\beta} &= \frac{xz}{c} \\
        \cos{\gamma} &= \frac{xy}{b}

Additional documentation can be found at `lammps <http://lammps.sandia.gov/doc/Section_howto.html?highlight=prism#howto-12>`_.
        """
        cdef int dim = self.dimension
        cdef double[::1] boxlo = <double[:dim]>self._domain.boxlo
        cdef double[::1] boxhi = <double[:dim]>self._domain.boxhi
        cdef double lx = boxhi[0] - boxlo[0]
        cdef double ly = boxhi[1] - boxlo[1]
        cdef double lz = boxhi[2] - boxlo[2]
        cdef double xy = self._domain.xy
        cdef double xz = self._domain.xz
        cdef double yz = self._domain.yz
        cdef double a = lx
        cdef double b = sqrt(ly**2 + xy**2)
        cdef double c = sqrt(lz**2 + xz**2 + yz**2)
        cdef double alpha = acos((xy*xz + ly*yz) / (b*c))
        cdef double beta = acos(xz / c)
        cdef double gamma = acos(xy / b)
        return (a, b, c), (alpha, beta, gamma)

    @property
    def lengths(self):
        """ See :py:func:`lengths_angles`

        :getter: return lengths (a, b, c)
        """
        return self.lengths_angles[0]

    @property
    def angles(self):
        """ See :py:func:`lengths_angles`

        :getter: returns angles (alpha, beta, gamma)
        """
        return self.lengths_angles[1]

    @property
    def volume(self):
        """ Volume of box given in lammps units 

        :getter: Returns volume of box given in lammps units 
        """
        cdef double vol = self._domain.xprd * self._domain.yprd
        if self.dimension == 2:
            return vol
        else: # dimension == 3
            return vol * self._domain.zprd
