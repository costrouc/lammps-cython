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

# Imports C-level symbols
cimport numpy as np
from libc.stdlib cimport malloc, free
from libc.math cimport sqrt, acos, cos

# Import Python-level symbols
from mpi4py import MPI
import numpy as np
from math import pi

cdef class Lammps:
    """LAMMPS base class represents the entire library.

..  py:function:: __init__(self, units='lj', style='atomic', comm=None, args=None)
    
    Initialize a Lammps object. 

    :param str units: units to use for simulation
    :param str style: atomic style to use for simulation
    :param comm: mpi4py comm object default value is MPI_COMM_WORLD
    :param listargs: list of command line args that would be supplied to normal lammps executable

    Possible values for the following arguments: 
     * args: `command-line arguments <http://lammps.sandia.gov/doc/Section_start.html#command-line-options>`_ 
     * units: `units <http://lammps.sandia.gov/doc/units.html>`_
     * style: `atom_style <http://lammps.sandia.gov/doc/atom_style.html>`_
    """
    cdef LAMMPS *_lammps
    cdef mpi.MPI_Comm _comm
    cdef public Box box
    cdef public System system
    cdef public Thermo thermo

    available_units = ['real', 'metal', 'si', 'cgs', 'electron', 'micro', 'nano', 'lj']
    
    def __cinit__(self, units='lj', style='atomic', comm=None, args=None):
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
        
        if args == None:
            args = ['python']
        else:
            args = ['python'] + args
            
        cdef int argc = len(args)
        cdef char** argv = <char**>args_to_cargv(args)

        self._lammps = new LAMMPS(argc, argv, self._comm)
        # don't free char** becuase used by lammps (they don't copy their strings!)

        # Set the units
        if units in Lammps.available_units:
            self.command("units {}".format(units))
        else:
            raise ValueError('units {} unknown units'.format(units))

        self.box = Box(self)
        self.system = System(self, style=style)
        self.thermo = Thermo(self)

    def __dealloc__(self):
        del self._lammps

    property __version__:
        """ Prints the version of LAMMPS 

        Format is <day><month><year> e.g. 7Dec15
        """
        def __get__(self):
            return (self._lammps.universe.version).decode('utf-8')

    def __repr__(self):
        rep = "<Lammps Style:{} Atoms:{:.2g} Lattice:[{:.1f}, {:.1f}, {:.1f}]>" 
        return rep.format(self.system.style, self.system.total, *self.box.lengths)

    def command(self, cmd):
        """Run a LAMMPS command

        :param str command: command for lammps to execute

        See lammps documentation for available `commands
        <http://lammps.sandia.gov/doc/Section_commands.html#individual-commands>`_.
        """
        self._lammps.input.one(cmd.encode('utf-8'))

    def file(self, filename):
        """ Run a LAMMPS input file

        :param string filename: filename of input file for LAMMPS

        This is equivalent to setting the -i command line flag.

..      todo:: learn how file() behaves with multiple invocations
        """
        self._lammps.input.file(filename.encode('utf-8'))

    def run(self, long steps):
        """ Runs the lammps simulation for N steps

        :param int steps: number of steps to run simulation equivalent to "run <steps>" command

        See lammps documentation for description of `run command
        <http://lammps.sandia.gov/doc/run.html>`_
        """
        self.command('run {}'.format(steps))

    def reset(self):
        """ Resets the lammps simulation

        Deletes all atoms, restores all settings to their default
        values, and frees all memory allocated by LAMMPS. Equivalent
        to the LAMMPS `clear command
        <http://lammps.sandia.gov/doc/clear.html>`_.
        """
        self.command('clear')

    property units:
        """ Units used in lammps simulation 

        :getter: Returns unit style used
        :type: string

        See `units <http://lammps.sandia.gov/doc/units.html>`_ for
        more information.
        """
        def __get__(self):
            
            return (self._lammps.update.unit_style).decode('utf-8')

    property dt:
        """ timestep size for run step in simulation **time units**
        
        :getter: Returns the timestep size
        :setter: Sets the timestep size
        :type: float
        """
        def __get__(self):
            return self._lammps.update.dt

        def __set__(self, double value):
            self._lammps.update.dt = value

    property time_step:
        """ current number of timesteps that have been run
        
        :getter: Returns the timestep number 
        :setter: Sets the timestep number
        :type: int
        """
        def __get__(self):
            return self._lammps.update.ntimestep

        def __set__(self, bigint value):
            self._lammps.update.reset_timestep(value)
    
    property time:
        """ total time that has elapsed from lammps runs in lammps **time units**
    
        :getter: Returns the total time
        :type: float
        """
        def __get__(self):
            return self._lammps.update.atime


cdef class Thermo:
    """ Computes thermodynamic properties of a group of particles

    You must first define a compute before you can extract
    thermodynamics properties of the current time step. Three computes
    are always created, named “thermo_temp”, “thermo_press”, and
    “thermo_pe” these are initialized in the output.cpp in
    LAMMPS. These are computes are made easy to access. All computes
    can be accessed via the comptutes dictionary.

..  py:function:: __init__(self, Lammps)
    
    Initialize a Thermo object.

    :param :py:class:`Lammps` lammps: Lammps object

..  py:attribute:: computes
    
    A dictionary of Computes. {id: :py:class:`Compute`}.
    """
    cdef Lammps lammps
    cdef public Compute temperature
    cdef public Compute pressure
    cdef public Compute potential_energy
    cdef public dict computes

    def __cinit__(self, Lammps lammps):
        """ Docstring in Thermo base class (sphinx can find doc when compiled) """
        self.lammps = lammps
        self.computes = dict()

        # Add computes automatically added by
        # Lammps (output.h)
        self.temperature = Compute(self.lammps, 'thermo_temp')
        self.computes.update({'thermo_temp': self.temperature})
        self.pressure = Compute(self.lammps, 'thermo_press')
        self.computes.update({'thermo_press': self.pressure})
        self.potential_energy = Compute(self.lammps, 'thermo_pe')
        self.computes.update({'thermo_pe': self.potential_energy})

    def add(self, id, style, group='all', args=None):
        """ Add a compute to LAMMPS

        :param str id: name of new lammps compute cannot conflict with existing compute ids
        :param str style: name of compute to add
        :param str group: name of compute group
        :param list args: additional args to supply to compute

        Equivalent lammps command: ``compute ID group-ID style args``

        See `compute <http://lammps.sandia.gov/doc/compute.html>`_ for
        more information on creating computes.
        """
        self.computes.update({id: Compute(self.lammps, id, style, group, args)})


cdef class Compute:
    """ Extract compute object property from LAMMPS system

    See `compute <http://lammps.sandia.gov/doc/compute.html>`_ for
    more information on available computes.

..  py:function:: __init__(self, Lammps)

    Initialize a Compute object.

    :param lammps: Lammps object
    """
    cdef Lammps lammps
    cdef COMPUTE* _compute

    def __cinit__(self, Lammps lammps, id, style=None, group='all', args=None):
        """Docstring in Compute base class (sphinx can find doc when compiled)"""
        self.lammps = lammps
        id = id.encode('utf-8')

        cdef int index = self.lammps._lammps.modify.find_compute(id)

        if index == -1:
            if style == None:
                raise ValueError("New computes require a style")

            cmd = (
                "compute {} {} {}"
            ).format(id, group, style)
            
            for arg in args:
                cmd += " {}".format(arg)

            self.lammps.command(cmd)
            index = self.lammps._lammps.modify.find_compute(id)

            if index == -1:
                raise Exception("Compute should have been created but wasn't")

        self._compute = self.lammps._lammps.modify.compute[index]

    def __repr__(self):
        rep = "<Compute:{} id:{}>"
        return rep.format(self.style, self.name)

    def modify(self, args):
        """ Modify a predefined LAMMPS compute

        :param list[str] args: args to modify lammps compute

        See `compute_modify
        <http://lammps.sandia.gov/doc/compute_modify.html>`_ for more
        information on args etc.

..      todo:: probably should just use lammps.command
        """
        cdef char** argv = args_to_cargv(args)
        cdef int argc = len(args)
        self._compute.modify_params(argc, argv)

    property name:
        """ Name of compute id
    
        :getter: Returns compute id
        :type: str  
        """
        def __get__(self):
            return (self._compute.id).decode("utf-8")

    property style:
        """ style of compute id 

        :getter: Returns compute style
        :type: str
        """
        def __get__(self):
            return (self._compute.style).decode("utf-8")

    property scalar:
        """ Scalar value of compute

        :getter: Returns scalar value of compute
        :type: float
        :raises NotImplementedError: no scalar function for this compute defined
        """
        def __get__(self):
            if self._compute.scalar_flag == 0:
                error_str = "Style {} does not have a scalar function"
                raise NotImplementedError(error_str.format(self.style))

            self._compute.compute_scalar()
            return self._compute.scalar

    property vector:
        """ vector value of compute
    
        :getter: Returns vector value of compute
        :type: numpy.ndarray
        :raises NotImplementedError: no vector function for this compute defined
        """
        def __get__(self):
            if self._compute.vector_flag == 0:
                error_str = "Style {} does not have a vector function"
                raise NotImplementedError(error_str.format(self.style))

            self._compute.compute_vector()

            cdef int N = self._compute.size_vector
            cdef double[::1] vector = <double[:N]>self._compute.vector
            return np.asarray(vector)


cdef class AtomType:
    """ Represents an atom type (denoted by int) in LAMMPS

    Provides a wrapper over atom types that is more user friendly

..  py:function:: __init__(self, Lammps lammps, int index)

    :param :py:class:`Lammps` lammps: lammps object 
    :param int index: the index of the atom

    Initialize an AtomType object.

..  py:attribute:: index
    
    (Read Only) Index of atom type within lammps

    :getter: Returns index of atom type
    :type: int
    """
    cdef Lammps lammps
    cdef readonly int index

    def __cinit__(self, Lammps lammps, int index):
        """ Docstring in AtomType base class (sphinx can find doc when compiled) """
        self.lammps = lammps

        if index < 1 or index > self.lammps._lammps.atom.ntypes:
            raise ValueError("Atom type index {} not in bounds".format(index))

        self.index = index

    property mass:
        """ Mass of the AtomType

        :getter: Returns the mass of AtomType
        :setter: Sets the mass of AtomType
        :type: float
        """
        def __get__(self):
            return self.lammps._lammps.atom.mass[self.index]

        def __set__(self, double value):
            self.lammps._lammps.atom.set_mass(self.index, value)


cdef class Atom:
    """ Represents a single atom in the LAMMPS simulation

    Since LAMMPS is a distributed system each atom is unique among all
    the processors.

..  py:function:: __init__(self, Lammps lammps, tag=None, lindex=None)

    :param :py:class:`Lammps` lammps: lammps object 
    :param int tag: tag number of particle 
    :param int lindex: the local index of the atom

    Initialize an System object. Must supply either a tag number or index

..  py:attribute:: tag

    Unique tag of atom within LAMMPS simulation

    :getter: Returns unqiue tag of atom
    :type: int
    """
    cdef Lammps lammps
    cdef long local_index
    cdef readonly long tag

    def __cinit__(self, Lammps lammps, tag=None, lindex=None):
        """ Docstring in System base class (sphinx can find doc when compiled) """
        self.lammps = lammps

        if lindex is not None:
            if lindex < 0 or lindex >= self.lammps.system.local:
                raise IndexError("index not withing local atoms index")
            self.local_index = lindex
            self.tag = self.lammps._lammps.atom.tag[self.local_index]
        else:
            raise NotImplementedError("currently need the local index")
            # self.tag = tag
            # self.index = self._atom.map(self.tag)

    property type:
        """ AtomType of local Atom

        :getter: Returns AtomType of the local atom
        :setter: sets the type of the local atom
        :type: AtomType or int
        """
        def __get__(self):
            cdef int itype = self.lammps._lammps.atom.type[self.local_index]
            return AtomType(self.lammps, itype)

        def __set__(self, value):
            if isinstance(value, AtomType):
                self.lammps._lammps.atom.type[self.local_index] = value.index
            elif isinstance(value, int):
                self.lammps._lammps.atom.type[self.local_index] = value
            else:
                raise TypeError("Setter must be AtomType or int")

    property position:
        """ Position of local Atom
        
        :getter: Returns the position of the local atom
        :type: np.ndarray[3]
        """
        def __get__(self):
            if self.lammps._lammps.atom.x == NULL:
                return None

            cdef double[::1] array = <double[:3]>&self.lammps._lammps.atom.x[0][self.local_index * 3]
            return np.asarray(array)

    property velocity:
        """ Velocity of local atom

        :getter: Returns the velocity of the local atom
        :type: np.ndarray[3]
        """
        def __get__(self):
            if self.lammps._lammps.atom.v == NULL:
                return None
            
            cdef double[::1] array = <double[:3]>&self.lammps._lammps.atom.v[0][self.local_index * 3]
            return np.asarray(array)

    property force:
        """ Force on local atom
        
        :getter: Returns the force of the local atom
        :type: np.ndarray[3]
        """
        def __get__(self):
            if self.lammps._lammps.atom.f == NULL:
                return None

            cdef double[::1] array = <double[:3]>&self.lammps._lammps.atom.f[0][self.local_index * 3]
            return np.asarray(array)

    property charge:
        """ Charge of local atom
        
        :getter: Returns the charge of the local atom
        :setter: Sets the charge of the local atom
        """
        def __get__(self):
            if self.lammps._lammps.atom.q == NULL:
                return None

            return self.lammps._lammps.atom.q[self.local_index]

        def __set__(self, double value):
            if self.lammps._lammps.atom.q == NULL:
                return

            self.lammps._lammps.atom.q[self.local_index] = value


cdef class System:
    """ Represents all the atoms in the LAMMPS simulation

    Since LAMMPS is a distributed system each processor has a local
    view of its Atoms.

..  py:function:: __init__(self, Lammps, style='atomic')
    
    Initialize a System object.

    :param :py:class:`Lammps` lammps: Lammps object
    :param str style: one the `atom_styles <http://lammps.sandia.gov/doc/atom_style.html>`_
    """
    cdef ATOM* _atom
    cdef Lammps lammps
    cdef unsigned int local_index # used by iter

    # Available styles (but not all properties are accessible)
    # More will be added in the future
    available_styles = [
        'angle', 'atomic', 'body', 'bond', 'charge'
        'dipole', 'electron', 'ellipsoid', 'full',
        'line', 'meso', 'molecular', 'peri', 'smd',
        'sphere', 'template', 'tri', 'wavepacket'
    ]

    def __cinit__(self, Lammps lammps, style='atomic'):
        """ Docstring in System base class (sphinx can find doc when compiled) """
        self.lammps = lammps
        self.local_index = 0

        # Set the atomic style
        if style in System.available_styles:
            self.lammps.command("atom_style {}".format(style))
        else:
            raise ValueError('style {} is an invalid style'.format(style))

    property style:
        """ The LAMMPS `atom_style <http://lammps.sandia.gov/doc/atom_style.html>`_ used in simulation

        :getter: Returns the `atom_style <http://lammps.sandia.gov/doc/atom_style.html>`_ used in simulation
        :type: str
        """
        def __get__(self):
            return (self.lammps._lammps.atom.atom_style).decode('utf-8')

    property total:
        """ Total number of atoms in LAMMPS simulation

        :getter: Returns the total number of atoms in LAMMPS simulation
        :type: int
        """
        def __get__(self):
            return self.lammps._lammps.atom.natoms

    property local:
        """ Local number of atoms stored on processor.

        :getter: Returns the local number of atoms specific to core
        :type: int

        An important concept in LAMMPS parallal calculations is that
        the atoms are distributed between all the processors.
        """
        def __get__(self):
            return self.lammps._lammps.atom.nlocal

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

        atom = Atom(self.lammps, lindex=self.local_index)
        self.local_index += 1
        return atom

    property tags:
        """ Tags associated with local atoms stored on processor

        :getter: Returns the local tags of atoms specific to core
        :type: nd.array[int](local)
        """
        def __get__(self):
            if self.lammps._lammps.atom.x == NULL:
                return None
            
            cdef size_t N = self.local
            cdef tagint[::1] array = <tagint[:N]>self.lammps._lammps.atom.tag
            return np.asarray(array)

    property types:
        """ Types associated with local atoms stored on processor

        :getter: Returns the local int types of atoms specific to core
        :type: nd.array[int](local)
        """
        def __get__(self):
            if self.lammps._lammps.atom.x == NULL:
                return None
            
            cdef size_t N = self.local
            cdef int[::1] array = <int[:N]>self.lammps._lammps.atom.type
            return np.asarray(array)

    property atom_types:
        """ Atom Types that are currently defined in LAMMPS simulation.

        :getter: Returns AtomTypes of system
        :type: list[AtomTypes]

..      note::
        
        Atom Types in LAMMPS start at an index of 1.
        """
        def __get__(self):
            atomtypes = []
            # AtomType int begins at 1 not 0.
            for i in range(1, self.lammps._lammps.atom.ntypes + 1):
                atomtypes.append(AtomType(self.lammps, i))
            return atomtypes

    property positions:
        """ Positions associated with local atoms stored on processor

        :getter: Returns the local positions of atoms specific to processor
        :type: np.array[double](local, 3)
        """
        def __get__(self):
            if self.lammps._lammps.atom.x == NULL:
                return None

            cdef size_t N = self.local
            cdef double[:, ::1] array = <double[:N, :3]>self.lammps._lammps.atom.x[0]
            return np.asarray(array)

    property velocities:
        """ Velocities associated with local atoms stored on processor

        :getter: Returns the local velocities of atoms specific to processor
        :type: np.array[double](local, 3)
        """
        def __get__(self):
            if self.lammps._lammps.atom.v == NULL:
                return None
            
            cdef size_t N = self.local
            cdef double[:, ::1] array = <double[:N, :3]>self.lammps._lammps.atom.v[0]
            return np.asarray(array)

    property forces:
        """ Forces associated with local atoms stored on processor

        :getter: Returns the local forces of atoms specific to processor
        :type: np.array[double](local, 3)
        """
        def __get__(self):
            if self.lammps._lammps.atom.f == NULL:
                return None
            
            cdef size_t N = self.local
            cdef double[:, ::1] arr = <double[:N, :3]>self.lammps._lammps.atom.f[0]
            return np.asarray(arr)

    property charges:
        """ Charges associated with local atoms stored on processor

        :getter: Returns the local charges of atoms specific to processor
        :type: np.array[double](local, 3)
        """
        def __get__(self):
            if self.lammps._lammps.atom.q == NULL:
                return None

            cdef size_t N = self.local
            cdef double[::1] vector = <double[:N]>self.lammps._lammps.atom.q
            return np.asarray(vector)

    def add(self, atom_type, position, remap=False, units='box'):
        """Create atom in LAMMPS system

        :param int atom_type: atom type id
        :param np.array[3] position: position of atom
        :param boolean remap: whether to remap atoms in box or not
        :param str units: units to use for positions

        essentially runs command: ``create_atoms <atom_type> single <x> <y> <z>``

        The lammps internals are hard to use. Yes there are some
        expensive mpi calls becuase of this (so this is slow for large
        systems). HACK

        See lammps `create_atoms
        <http://lammps.sandia.gov/doc/create_atoms>`_. for documentation.
        """
        if len(position) != 3:
            raise ValueError('position array must be 3 values')

        cmd = (
            'create_atoms {!s} single {!s} {!s} {!s} units {!s} remap {!s}'
        ).format(atom_type, position[0], position[1], position[2], units,
                 'yes' if remap else 'no')
        self.lammps.command(cmd)


cdef class Box:
    """ Represents the shape of the simulation cell.

..  py:function:: __init__(self, :py:class:`Lammps` lammps)

    Initialize a Box object.

    :param :py:class:`Lammps` lammps: Lammps object
    """
    cdef Lammps lammps

    def __cinit__(self, Lammps lammps):
        """ Docstring in Box base class (sphinx can find doc when compiled) """
        self.lammps = lammps

    def from_lattice_const(self, atom_types, lengths, angles=None, region_id='lammps-box'):
        """ Create Unit cell (can only be run once)

        :param int atom_types: number of atom types for simulation (I know weird place to require)
        :param list lengths: lattice constants of unit cell
        :param list angles: angles of unit cell in radians
        :param str region_id: name of region for lammps simulation

        Essentially runs:
        region region-id prism xlo xhi ylo yhi zlo zhi xy xz yz
        create_box N region-ID keyword value ...

..      todo::
        This functions needs more thought. It does too much but within
        the LAMMPS C++ function calls I can't see a nice way to do it.
        """
        if angles == None:
            angles = [pi/2., pi/2., pi/2.]

        (lx, ly, lz), (xy, xz, yz) = lattice_const_to_lammps_box(lengths, angles)
        self.lammps.command((
            "region {} prism 0.0 {} 0.0 {} 0.0 {} {} {} {}"
        ).format(region_id, lx, ly, lz, xy, xz, yz))

        self.lammps.command((
            "create_box {} {}"
        ).format(atom_types, region_id))

    property dimension:
        """ The dimension of the lammps run either 2D or 3D

        :getter: Returns the dimension of lammps simulation (2 or 3)
        :type: int
        """
        def __get__(self):
            return self.lammps._lammps.domain.dimension

    property lohi:
        """ LAMMPS box description of boxhi and boxlo

        :return: dictionary of lower and upper in each dimension
        :rtype: dict

        Additional documentation can be found at `lammps
        <http://lammps.sandia.gov/doc/Section_howto.html?highlight=prism#howto-12>`_. For
        example a return dictionary would be

..      code-block:: python

        lohi = {
            'boxlo': np.array([0.0, 0.0, 0.0]),
            'boxhi': np.array([10.0, 10.0, 10.0])
        }
        """
        def __get__(self):
            cdef int dim = self.dimension
            cdef double[::1] boxlo = <double[:dim]>self.lammps._lammps.domain.boxlo
            cdef double[::1] boxhi = <double[:dim]>self.lammps._lammps.domain.boxhi
            return {
                'boxlo': np.array(boxlo),
                'boxhi': np.array(boxhi)
            }

    property tilts:
        """ LAMMPS box description of xy xz yz

        :return: dictionary of xy xz yz
        :rtype: dict

        Additional documentation can be found at `lammps
        <http://lammps.sandia.gov/doc/Section_howto.html?highlight=prism#howto-12>`_.
        For example a return dictionary would be

..      code-block:: python

        tilts = {
           'xy': 0.0,
           'xz': 0.0,
           'yz': 0.0
        }

..      todo::
        
        how to handle 2d? Needs to be addressed throughout
        """
        def __get__(self):
            cdef double xy = self.lammps._lammps.domain.xy
            cdef double xz = self.lammps._lammps.domain.xz
            cdef double yz = self.lammps._lammps.domain.yz
            return {
                'xy': xy, 
                'xz': xz, 
                'yz': yz
            }

    property lengths_angles:
        """ Calculates the lengths and angles from lammps box

        :getter: return lengths (a, b, c) and angles (alpha, beta, gamma)
        :type: tuple[double](3), tuple[double](3)

        Additional documentation can be found at `lammps
        <http://lammps.sandia.gov/doc/Section_howto.html?highlight=prism#howto-12>`_. See
        implementation :py:func:`lammps_box_to_lattice_const`.
        """
        def __get__(self):
            cdef int dim = self.dimension
            cdef double[::1] boxlo = <double[:dim]>self.lammps._lammps.domain.boxlo
            cdef double[::1] boxhi = <double[:dim]>self.lammps._lammps.domain.boxhi
            cdef double lx = boxhi[0] - boxlo[0]
            cdef double ly = boxhi[1] - boxlo[1]
            cdef double lz = boxhi[2] - boxlo[2]
            cdef double xy = self.lammps._lammps.domain.xy
            cdef double xz = self.lammps._lammps.domain.xz
            cdef double yz = self.lammps._lammps.domain.yz
            return lammps_box_to_lattice_const((lx, ly, lz), (xy, xz, yz))

    property lengths:
        """ See :py:func:`lengths_angles`

        :getter: return lengths (a, b, c)
        :type: tuple[float](3)
        """
        def __get__(self):
            return self.lengths_angles[0]

    property angles:
        """ See :py:func:`lengths_angles`

        :getter: returns angles in radians (alpha, beta, gamma)
        :type: tuple[float](3)
        """
        def __get__(self):
            return self.lengths_angles[1]

    property volume:
        """ Volume of box given in lammps units 

        :getter: Returns volume of box given in lammps units 
        :type: float
        """
        def __get__(self):
            cdef double xprd = self.lammps._lammps.domain.xprd
            cdef double yprd = self.lammps._lammps.domain.yprd
            cdef double vol = xprd * yprd
            if self.dimension == 2:
                return vol
            else:
                return vol * self.lammps._lammps.domain.zprd

def lattice_const_to_lammps_box(lengths, angles):
    """ Converts lattice constants to lammps box coorinates(angles in radians) 

..  math::

    lx &= boxhi[0] - boxlo[0] = a \\
    ly &= boxhi[1] - boxlo[1] = \sqrt{b^2 - xy^2} \\
    lz &= boxhi[2] - boxlo[2] = \sqrt{c^2 - xz^2 - yz^2}
    xy &= b \cos{\gamma} \\
    xz &= c \cos{\beta}  \\
    yz &= \frac{b c \cos{\alpha} - xy xz}{ly}
    """
    a, b, c = lengths
    alpha, beta, gamma = angles
    cdef lx = a
    cdef xy = b * cos(gamma)
    cdef xz = c * cos(beta)
    cdef ly = sqrt(b**2 - xy**2)
    cdef yz = (b*c*cos(alpha) - xy*xz)/ly
    cdef lz = sqrt(c**2 - xz**2 - yz**2)
    return (lx, ly, lz), (xy, xz, yz)


def lammps_box_to_lattice_const(lengths, tilts):
    """ Converts lammps box coordinates to lattice constants

..  math::

    a &= lx \\
    b &= \sqrt{ly^2 + xy^2} \\
    c &= \sqrt{lz^2 + xz^2 + yz^2} \\
    \cos{\alpha} &= \frac{xy xz + ly yz}{b c} \\
    \cos{\beta} &= \frac{xz}{c} \\
    \cos{\gamma} &= \frac{xy}{b}
    """
    lx, ly, lz = lengths
    xy, xz, yz = tilts
    cdef double a = lx
    cdef double b = sqrt(ly**2 + xy**2)
    cdef double c = sqrt(lz**2 + xz**2 + yz**2)
    cdef double alpha = acos((xy*xz + ly*yz) / (b*c))
    cdef double beta = acos(xz / c)
    cdef double gamma = acos(xy / b)
    return (a, b, c), (alpha, beta, gamma)


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
