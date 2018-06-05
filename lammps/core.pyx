#!python
# cython: embedsignature=True

"""LAMMPS Python interface (the GOOD parts)

This interface is inspired by HOOMD and tries its best to look
similar. I try to include only orthogonal functions (e.g. make
there only be ONE way to do something).

.. todo:: Features
 - Group support
 - integrator support
 - lammps units support (maybe integrate with units package)
"""

include "core.pyd"

# Imports C-level symbols
cimport numpy as np
from libc.stdlib cimport malloc, free
from libc.string cimport strcpy, strlen
from libc.math cimport sqrt, acos, cos
from libcpp cimport bool

# Import Python-level symbols
from mpi4py import MPI
import numpy as np
from math import pi
import warnings


class LammpsDataTypes:
    INT = np.intc
    SMALL_INT = np.intc
    IMAGE_INT = np.intc
    TAG_INT = np.intc
    BIG_INT = np.int64
    DOUBLE = np.float64


cdef class LammpsError(Exception):
    def __cinit__(self, message):
        self.message = message


cdef class LammpsNormalError(LammpsError):
    pass


cdef class LammpsAbortError(LammpsError):
    pass


cdef class Lammps:
    """LAMMPS base class represents the entire library.

    Parameters
    ----------
    units : :obj:`str`
        units to use for simulation default lj
    style : :obj:`str`
        atomic style to use for simulation. default atomic
    comm : :obj:`mpi4py.Comm`
        default value is MPI_COMM_WORLD
    args : :obj:`list[str]`
        command line args that would be supplied to normal lammps executable

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

    AVAILABLE_UNITS = {'real', 'metal', 'si', 'cgs', 'electron', 'micro', 'nano', 'lj'}

    def __cinit__(self, str units='lj', str style='atomic', comm=None, args=None):
        """Docstring in Lammps base class (sphinx can find doc when compiled)"""
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
        if units in Lammps.AVAILABLE_UNITS:
            self.command("units {}".format(units))
        else:
            raise ValueError('units {} unknown units'.format(units))

        self.box = Box(self)
        self.system = System(self, style=style)
        self.thermo = Thermo(self)

        # TODO hack run file if -i flag provided.
        # Eventually abstract away command line arguments
        # main.cpp calls input->file() automatically
        if '-in' in args or '-i' in args:
            self._lammps.input.file()

    def __dealloc__(self):
        del self._lammps

    def check_error(self):
        """Check for error from lammps"""
        cdef char* error_message = self._lammps.error.get_last_error()
        if error_message == NULL:
            return None
        cdef ErrorType error_type = self._lammps.error.get_last_error_type()
        cdef str error_message_utf = error_message.decode('utf-8')
        self._lammps.error.set_last_error(NULL, ErrorType.ERROR_NONE)
        if error_type == ErrorType.ERROR_NORMAL:
            raise LammpsNormalError(error_message_utf)
        elif error_type == ErrorType.ERROR_ABORT:
            raise LammpsAbortError(error_message_utf)

    @property
    def __version__(self):
        """ Prints the version of LAMMPS

        Format is "<day> <month> <year>" e.g. 16 Mar 2018
        """
        return self._lammps.universe.version.decode('utf-8')

    def __repr__(self):
        rep = "<Lammps Style:{} Atoms:{:.2g} Lattice:[{:.1f}, {:.1f}, {:.1f}]>"
        return rep.format(self.system.style, self.system.total, *self.box.lengths)

    def command(self, str cmd):
        """Run a LAMMPS command

        Parameters
        ----------
        command : :obj:`str`
            command for lammps to execute

        Returns
        -------
        str
            output from running command

        See lammps documentation for available `commands
        <http://lammps.sandia.gov/doc/Section_commands.html#individual-commands>`_.
        """
        cdef char* result = lammps_command(self._lammps, cmd.encode('utf-8'))
        self.check_error()
        if result == NULL:
            return None
        return result.decode('utf-8')

    def file(self, str filename):
        """Read in a file as sequential list of commands

        Parameters
        ----------
        filename : :obj:`str`
            filename of input file for LAMMPS

        Can be invoked more than once
        """
        lammps_file(self._lammps, filename.encode('utf-8'))
        self.check_error()

    def run(self, long steps, bool pre=True, bool post=True):
        """ Runs the lammps simulation for N steps

        Parameters
        ----------
        steps : :obj:`int`
             number of steps to run simulation equivalent to "run <steps>" command
        pre : :obj:`bool`
             before run create neighbor lists, compute forces, and imposes fix contraints. default True
        post : :obj:`bool`
             print full timing summary

        See lammps documentation for description of `run command
        <http://lammps.sandia.gov/doc/run.html>`_
        """
        self.command('run {} pre {} post {}'.format(steps,
                                     "yes" if pre  else "no",
                                     "yes" if post else "no"))

    def reset(self):
        """Resets (clear) the lammps simulation

        Deletes all atoms, restores all settings to their default
        values, and frees all memory allocated by LAMMPS. Equivalent
        to the LAMMPS `clear command
        <http://lammps.sandia.gov/doc/clear.html>`_. Some settings are
        not affected. These include: working directory, log file
        status, echo status, and input script variables.
        """
        self.command('clear')

    @property
    def units(self):
        """ Unit style used in lammps simulation

        See `units <http://lammps.sandia.gov/doc/units.html>`_ for
        more information.
        """
        return self._lammps.update.unit_style.decode('utf-8')

    @property
    def dt(self):
        """ Timestep size for run step in simulation **time units**"""
        return self._lammps.update.dt

    @dt.setter
    def dt(self, double value):
        self._lammps.update.dt = value

    @property
    def time_step(self):
        """ current number of timesteps that have been run"""
        return self._lammps.update.ntimestep

    @time_step.setter
    def time_step(self, bigint value):
        self._lammps.update.reset_timestep(value)

    @property
    def time(self):
        """total time that has elapsed from lammps runs **time units**"""
        return self._lammps.update.atime


cdef class Thermo:
    """ Computes thermodynamic properties of a group of particles

    You must first define a compute before you can extract
    thermodynamics properties of the current time step. Three computes
    are always created, named “thermo_temp”, “thermo_press”, and
    “thermo_pe” these are initialized in the output.cpp in
    LAMMPS. These computes are made easy to access. All computes can
    be accessed via the comptutes dictionary.

    Parameters
    ----------
    lammps : :class:`Lammps`
        lammps object
    """
    cdef Lammps lammps
    cdef public Compute temperature
    cdef public Compute pressure
    cdef public Compute potential_energy
    cdef public dict _computes

    def __cinit__(self, Lammps lammps):
        """ Docstring in Thermo base class (sphinx can find doc when compiled) """
        self.lammps = lammps
        self._computes = dict()

        # Add computes automatically added by
        # Lammps (output.h)
        self.temperature = Compute(self.lammps, 'thermo_temp')
        self._computes.update({'thermo_temp': self.temperature})
        self.pressure = Compute(self.lammps, 'thermo_press')
        self._computes.update({'thermo_press': self.pressure})
        self.potential_energy = Compute(self.lammps, 'thermo_pe')
        self._computes.update({'thermo_pe': self.potential_energy})

    def add(self, str id, str style, str group='all', args=None):
        """ Add a compute to LAMMPS

        Parameters
        ----------
        id : str
           name of new LAMMPS compute must be unique
        style : str
           name of compute to add
        group : str
           name of compute group. default 'all'.
        args : list[str]
           additional args to supply to compute

        Equivalent lammps command: ``compute ID group-ID style args``

        See `compute <http://lammps.sandia.gov/doc/compute.html>`_ for
        more information on creating computes.
        """
        self._computes.update({id: Compute(self.lammps, id, style, group, args)})

    @property
    def computes(self):
        """Dictionary of available computes

        """
        return self._computes


cdef class Compute:
    """ Extract compute object property from LAMMPS system

    See `compute <http://lammps.sandia.gov/doc/compute.html>`_ for
    more information on available computes. This is how properties of
    systems can be exctracted without I/O. See list of `available
    computes
    <http://lammps.sandia.gov/doc/Section_commands.html#cmd-5>`

    Parameters
    ----------
    lammps : :class:`lammps.Lammps`
          Lammps object
    id : str
          name id for lammps compute must be unique
    style : :obj:`str`
          compute style to use see `LAMMPS documentation Compute
          styles
          <http://lammps.sandia.gov/doc/Section_commands.html#cmd-5>`_
    group : :obj:`str`
          group to apply compute on. default "all" must create group
          if using other group. See `group command
          <http://lammps.sandia.gov/doc/group.html>`_ if creating own
          group
    args : list[str]
         list of additional arguments to supply.
    """
    cdef Lammps lammps
    cdef COMPUTE* _compute

    def __cinit__(self, Lammps lammps, str id, style=None, group='all', args=None):
        """Docstring in Compute base class (sphinx can find doc when compiled)"""
        self.lammps = lammps

        cdef bytes bytes_id = id.encode('utf-8')

        if args is None:
            args = []

        cdef int index = self.lammps._lammps.modify.find_compute(bytes_id)
        if index == -1: # Compute Id does not currently exist
            if style == None:
                raise ValueError("New computes require a style")

            cmd = (
                "compute {} {} {}"
            ).format(id, group, style)

            for arg in args:
                cmd += " {}".format(arg)

            self.lammps.command(cmd)
            index = self.lammps._lammps.modify.find_compute(bytes_id)

            if index == -1:
                # TODO: Make Lammps exceptions and errors
                raise LammpsNormalError("Compute should have been created but wasn't")

        self._compute = self.lammps._lammps.modify.compute[index]

    def __repr__(self):
        rep = "<Compute:{} id:{}>"
        return rep.format(self.style, self.name)

    def modify(self, args):
        """ Modify a predefined LAMMPS compute

        Parameters
        ----------
        args : list[str]
             args to modify lammps compute

        See `compute_modify
        <http://lammps.sandia.gov/doc/compute_modify.html>`_ for more
        information on args etc.
        """
        cdef char** argv = args_to_cargv(args)
        cdef int argc = len(args)
        self._compute.modify_params(argc, argv)

    @property
    def name(self):
        """Name of compute id"""
        return (self._compute.id).decode("utf-8")

    @property
    def style(self):
        """ style of compute id"""
        return (self._compute.style).decode("utf-8")

    @property
    def scalar(self):
        """ Scalar value of compute

        Returns
        -------
        float
              scalar value of compute

        Raises
        ------
        NotImplementedError
              no scalar function for this compute defined
        """
        if self._compute.scalar_flag == 0:
            error_str = "Style {} does not have a scalar function"
            raise NotImplementedError(error_str.format(self.style))

        self._compute.compute_scalar()
        return self._compute.scalar

    @property
    def vector(self):
        """ vector value of compute

        np.ndarray
              vector value of compute

        Raises
        ------
        NotImplementedError
              no vector function for this compute defined
        """
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

    Parameters
    ----------
    lammps : :class:`lammps.Lammps`
        Lammps object
    index : int
        the index of the atom
    """
    cdef Lammps lammps
    cdef readonly int _index

    def __cinit__(self, Lammps lammps, int index):
        """ Docstring in AtomType base class (sphinx can find doc when compiled) """
        self.lammps = lammps

        if index < 1 or index > self.lammps._lammps.atom.ntypes:
            raise ValueError("Atom type index {} not in bounds".format(index))

        self._index = index

    @property
    def index(self):
        """Index of atom type"""
        return self._index

    @property
    def mass(self):
        """Mass of the AtomType"""
        return self.lammps._lammps.atom.mass[self.index]

    @mass.setter
    def mass(self, double value):
        self.lammps._lammps.atom.set_mass("python", 1, self.index, value)


cdef class Atom:
    """ Represents a single atom in the LAMMPS simulation

    Since LAMMPS is a distributed system each atom is unique among all
    the processors.

    Parameters
    ----------
    lammps : :class:`lammps.core.Lammps`
         Lammps object
    tag : int
         tag id number of particle
    lindex : int
         lindex: the local index of the atom
    """
    cdef Lammps lammps
    cdef long local_index
    cdef readonly long _tag

    def __cinit__(self, Lammps lammps, tag=None, lindex=None):
        """ Docstring in System base class (sphinx can find doc when compiled) """
        self.lammps = lammps

        if lindex is not None:
            if lindex < 0 or lindex >= self.lammps.system.local:
                raise IndexError("index not withing local atoms index")
            self.local_index = lindex
            self._tag = self.lammps._lammps.atom.tag[self.local_index]
        else:
            raise NotImplementedError("currently need the local index")
            # self.tag = tag
            # self.index = self._atom.map(self.tag)

    @property
    def tag(self):
        """Unique tag of atom within LAMMPS simulation"""
        return self._tag

    @property
    def type(self):
        """AtomType of local Atom"""
        cdef int itype = self.lammps._lammps.atom.type[self.local_index]
        return AtomType(self.lammps, itype)

    @type.setter
    def type(self, value):
        if isinstance(value, AtomType):
            self.lammps._lammps.atom.type[self.local_index] = value.index
        elif isinstance(value, int):
            self.lammps._lammps.atom.type[self.local_index] = value
        else:
            raise TypeError("Setter must be AtomType or int")

    @property
    def position(self):
        """ Position of local Atom"""
        if self.lammps._lammps.atom.x == NULL:
            return None

        cdef double[::1] array = <double[:3]>&self.lammps._lammps.atom.x[0][self.local_index * 3]
        return np.asarray(array)

    @property
    def velocity(self):
        """ Velocity of local atom"""
        if self.lammps._lammps.atom.v == NULL:
            return None

        cdef double[::1] array = <double[:3]>&self.lammps._lammps.atom.v[0][self.local_index * 3]
        return np.asarray(array)

    @property
    def force(self):
        """ Force on local atom"""
        if self.lammps._lammps.atom.f == NULL:
            return None

        cdef double[::1] array = <double[:3]>&self.lammps._lammps.atom.f[0][self.local_index * 3]
        return np.asarray(array)

    @property
    def charge(self):
        """Charge of local atom"""
        if self.lammps._lammps.atom.q == NULL:
            return None

        return self.lammps._lammps.atom.q[self.local_index]

    @charge.setter
    def charge(self, double value):
        if self.lammps._lammps.atom.q == NULL:
            return

        self.lammps._lammps.atom.q[self.local_index] = value


cdef class System:
    """Represents all the atoms and respective properties in the LAMMPS simulation

    Since LAMMPS is a distributed system each processor has a local
    view of its Atoms. This class hopes to provide both and local and
    global view of the properties.

    Parameters
    ----------
    lammps : :class:`lammps.core.Lammps`
         LAMMPS Object
    style : :obj:`str`
         one the `atom_styles <http://lammps.sandia.gov/doc/atom_style.html>`_
    """
    cdef ATOM* _atom
    cdef Lammps lammps
    cdef unsigned int local_index # used by iter

    # Available styles (but not all properties are accessible)
    # More will be added in the future
    AVAILABLE_ATOM_STYLES = [
        'angle', 'atomic', 'body', 'bond', 'charge'
        'dipole', 'electron', 'ellipsoid', 'full',
        'line', 'meso', 'molecular', 'peri', 'smd',
        'sphere', 'template', 'tri', 'wavepacket'
    ]

    # from atom.cpp extract(name)
    ATOM_STYLE_PROPERTIES = {
        'mass': (LammpsDataTypes.DOUBLE, 1), # mass of particle (mass units)
        'id': (LammpsDataTypes.TAG_INT, 1), # integer ID of atom
        'type': (LammpsDataTypes.INT, 1), # type of atom (1-Ntype)
        'mask': (LammpsDataTypes.INT, 1),
        'image': (LammpsDataTypes.IMAGE_INT, 1),
        'x': (LammpsDataTypes.DOUBLE, 3), # atom position (position units)
        'v': (LammpsDataTypes.DOUBLE, 3), # atom velocity (velocity units)
        'f': (LammpsDataTypes.DOUBLE, 3), # atom force (force units)
        'molecule': (LammpsDataTypes.TAG_INT, 1), # integer ID of molecule the atom belongs to
        'q': (LammpsDataTypes.DOUBLE, 1), # charge on atom (charge units)
        'mu': (LammpsDataTypes.DOUBLE, 3), # x,y,z components of dipole moment of atom (dipole units)
        'omega': (LammpsDataTypes.DOUBLE, 3),
        'angmom': (LammpsDataTypes.DOUBLE, 3),
        'torque': (LammpsDataTypes.DOUBLE, 3),
        'radius': (LammpsDataTypes.DOUBLE, 1),
        'rmass': (LammpsDataTypes.DOUBLE, 1),
        'ellipsoid': (LammpsDataTypes.DOUBLE, 1),
        'line': (LammpsDataTypes.DOUBLE, 1),
        'tri': (LammpsDataTypes.DOUBLE, 1),
        # Peri Package
        'vfrac': (LammpsDataTypes.DOUBLE, 3),
        's0': (LammpsDataTypes.DOUBLE, 1),
        'x0': (LammpsDataTypes.DOUBLE, 3),
        # USER-EFF & USER-AWPMD
        'spin': (LammpsDataTypes.INT, 1),
        'eradius': (LammpsDataTypes.DOUBLE, 1), # electron radius (or fixed-core radius)
        'ervel': (LammpsDataTypes.DOUBLE, 1),
        'erforce': (LammpsDataTypes.DOUBLE, 1),
        'ervelforce': (LammpsDataTypes.DOUBLE, 1),
        'cs': (LammpsDataTypes.DOUBLE, 1),
        'csforce': (LammpsDataTypes.DOUBLE, 1),
        'vforce': (LammpsDataTypes.DOUBLE, 1),
        'etag': (LammpsDataTypes.INT),
        # USER-SPH
        'rho': (LammpsDataTypes.DOUBLE, 1), # density (need units) for SPH particles
        'drho': (LammpsDataTypes.DOUBLE, 1),
        'e': (LammpsDataTypes.DOUBLE, 1),
        'de': (LammpsDataTypes.DOUBLE, 1),
        'cv': (LammpsDataTypes.DOUBLE, 1),
        'vest': (LammpsDataTypes.DOUBLE, 3),
        # USER-SMD
        'contact_radius': (LammpsDataTypes.DOUBLE, 1),
        'smd_data_9': (LammpsDataTypes.DOUBLE, 3),
        'smd_stress': (LammpsDataTypes.DOUBLE, 3),
        'eff_plastic_strain': (LammpsDataTypes.DOUBLE, 1),
        'eff_plastic_strain_rate': (LammpsDataTypes.DOUBLE, 1),
        'damage': (LammpsDataTypes.DOUBLE, 1),
        # USER-DPD
        'dpdTheta': (LammpsDataTypes.DOUBLE, 1),
        # USER-MESO
        'edpd_temp': (LammpsDataTypes.DOUBLE, 1)
    }

    def __cinit__(self, Lammps lammps, style='atomic'):
        """ Docstring in System base class (sphinx can find doc when compiled) """
        self.lammps = lammps
        self.local_index = 0

        # Set the atomic style
        if style in System.AVAILABLE_ATOM_STYLES:
            self.lammps.command("atom_style {}".format(style))
        else:
            raise ValueError('style {} is an invalid style'.format(style))

    @property
    def style(self):
        """The LAMMPS `atom_style <http://lammps.sandia.gov/doc/atom_style.html>`_ used in simulation"""
        return self.lammps._lammps.atom.atom_style.decode('utf-8')

    @property
    def total(self):
        """Total number of atoms in LAMMPS simulation"""
        return self.lammps._lammps.atom.natoms

    @property
    def local_total(self):
        """ Local number of atoms stored on processor

        An important concept in LAMMPS parallal calculations is that
        the atoms are distributed between all the processors.
        """
        return self.lammps._lammps.atom.nlocal

    def __len__(self):
        return self.local_total

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

    @property
    def tags(self):
        """Tags associated with local atoms stored on processor"""
        return self.global_gather_property_ordered('id')

    @property
    def types(self):
        """Atom type ids associated with global atoms"""
        return self.global_gather_property_ordered('type')

    @property
    def atom_types(self):
        """AtomTypes in lammps system return list of AtomTypes"""
        atomtypes = []
        # AtomType int begins at 1 not 0.
        for i in range(1, self.lammps._lammps.atom.ntypes + 1):
            atomtypes.append(AtomType(self.lammps, i))
        return atomtypes

    @property
    def positions(self):
        """Positions associated with global atoms sorted by atom tag id

..      note::

        Reseting the atom positions is dangerous since it remaps the
        atoms to processors. I HIGHLY recommend create_atoms instead.
        """
        return self.global_gather_property_ordered('x')

    @positions.setter
    def positions(self, double[:, :] values):
        self.global_scatter_property_ordered('x', values)

    @property
    def velocities(self):
        """Velocities associated with global atoms sorted by atom tag id"""
        return self.global_gather_property_ordered('v')

    @velocities.setter
    def velocities(self, double[:, :] values):
        self.global_scatter_property_ordered('v', values)

    @property
    def forces(self):
        """Forces associated with global atoms sorted by atom tag id"""
        return self.global_gather_property_ordered('f')

    @forces.setter
    def forces(self, double[:, :] values):
        self.global_scatter_property_ordered('f', values)

    @property
    def charges(self):
        """Charges associated with global atoms sorted by tag id"""
        return self.global_gather_property_ordered('q')

    @charges.setter
    def charges(self, double[:, :] values):
        self.global_scatter_property_ordered('q', values)

    def global_gather_property_ordered(self, str name):
        """Gather globally system property to single processor. Sorted by atom id.

        Available properties are in :var:`System.ATOM_STYLE_PROPERTIES`
        """
        if name not in self.ATOM_STYLE_PROPERTIES:
            raise ValueError('atom system property %s does not exist' % name)
        atom_style_type, atom_style_count = self.ATOM_STYLE_PROPERTIES[name]
        data = np.zeros((self.total, atom_style_count), dtype=atom_style_type)
        if self.total == 0: # if no atoms just return empty array
            return data

        if atom_style_type == np.intc:
            self._global_gather_property_ordered_int(name.encode('utf-8'), atom_style_count, data)
        elif atom_style_type == np.float64:
            self._global_gather_property_ordered_double(name.encode('utf-8'), atom_style_count, data)
        else:
            raise TypeError('property %s type %s not recognized' % (name, atom_style_type))
        self.lammps.check_error()
        return data

    def _global_gather_property_ordered_int(self, char *name, int atom_style_count, int[:, :] data):
        lammps_gather_atoms(self.lammps._lammps, name, 0, atom_style_count, &data[0][0])

    def _global_gather_property_ordered_double(self, char *name, int atom_style_count, double[:, :] data):
        lammps_gather_atoms(self.lammps._lammps, name, 1, atom_style_count, &data[0][0])

    def global_gather_property_unordered(self, str name):
        """Gather globally system property to single processor. Data is not sorted.

        Available properties are in :var:`System.ATOM_STYLE_PROPERTIES`
        """
        if name not in self.ATOM_STYLE_PROPERTIES:
            raise ValueError('atom system property %s does not exist' % name)
        atom_style_type, atom_style_count = self.ATOM_STYLE_PROPERTIES[name]
        data = np.zeros((self.total, atom_style_count), dtype=atom_style_type)
        if self.total == 0: # if no atoms just return empty array
            return data

        if atom_style_type == np.intc:
            self._global_gather_property_unordered_int(name.encode('utf-8'), atom_style_count, data)
        elif atom_style_type == np.float64:
            self._global_gather_property_unordered_double(name.encode('utf-8'), atom_style_count, data)
        else:
            raise TypeError('property %s type %s not recognized' % (name, atom_style_type))
        self.lammps.check_error()
        return data

    def _global_gather_property_unordered_int(self, char *name, int atom_style_count, int[:, :] data):
        lammps_gather_atoms_concat(self.lammps._lammps, name, 0, atom_style_count, &data[0][0])

    def _global_gather_property_unordered_double(self, char *name, int atom_style_count, double[:, :] data):
        lammps_gather_atoms_concat(self.lammps._lammps, name, 1, atom_style_count, &data[0][0])

    def global_scatter_property_ordered(self, str name, data):
        """Scatter globally system property to all processors. Sorted by atom id.

        Available properties are in :var:`System.ATOM_STYLE_PROPERTIES`

        I HIGHLY recommend to not set atom positions with this
        method. It will result in lost atoms instead use create_atoms
        if not used properly
        - https://sourceforge.net/p/lammps/mailman/message/35842978/
        """
        if name not in self.ATOM_STYLE_PROPERTIES:
            raise ValueError('atom system property %s does not exist' % name)
        elif name == 'x':
            warnings.warn('setting atom positions using scatter may change processor ownership this can easily lead to lost atoms')

        atom_style_type, atom_style_count = self.ATOM_STYLE_PROPERTIES[name]
        if len(data) != self.total:
            raise ValueError('number of atoms must match data first dimmension')

        if self.lammps._lammps.atom.map_style == 0:
            raise ValueError('cannot set atom properties if map style is None set "atom_modify"')
        elif self.lammps._lammps.atom.tag_enable == 0:
            raise ValueError('tags must be enabled?')
        elif self.lammps._lammps.atom.tag_consecutive() == 0:
            raise ValueError('tags must be consecutive?')

        if atom_style_type == np.intc:
            self._global_scatter_property_ordered_int(name.encode('utf-8'), atom_style_count, data)
        elif atom_style_type == np.float64:
            self._global_scatter_property_ordered_double(name.encode('utf-8'), atom_style_count, data)
        else:
            raise TypeError('property %s type %s not recognized' % (name, atom_style_type))
        self.lammps.check_error()

    def _global_scatter_property_ordered_int(self, char *name, int atom_style_count, int[:, :] data):
        lammps_scatter_atoms(self.lammps._lammps, name, 0, atom_style_count, &data[0][0])

    def _global_scatter_property_ordered_double(self, char *name, int atom_style_count, double[:, :] data):
        lammps_scatter_atoms(self.lammps._lammps, name, 1, atom_style_count, &data[0][0])

    def create_atoms(self, int[:] atom_types, double[:, :] positions not None, double[:, :] velocities=None, start_index=None, bool wrap_atoms=True):
        """Create atoms for LAMMPS calculation.

        If atoms already exist on the system you will need to change
        the atom "start_id_index" to greater than the max value.

        Parameters
        ----------
        atom_types : list[int], np.array[int]
            vector of atom types
        positions : np.array[Nx3]
            array of atom positions
        velocities : np.array[Nx3]
            array of atom velocities. default set all velocities to 0
        start_index : int
            index to start atom ids
        """
        start_index = start_index or self.total + 1
        cdef int wrap_atoms_int = 1 if wrap_atoms else 0
        cdef int num_atoms = len(atom_types)
        cdef int[:] atom_ids = np.arange(start_index, len(atom_types) + start_index, dtype=np.intc)
        if velocities is None:
            lammps_create_atoms(self.lammps._lammps, num_atoms,
                                &atom_ids[0], &atom_types[0],
                                &positions[0, 0], NULL,
                                NULL, wrap_atoms_int)
        else:
            lammps_create_atoms(self.lammps._lammps, num_atoms,
                                &atom_ids[0], &atom_types[0],
                                &positions[0, 0], &velocities[0, 0],
                                NULL, wrap_atoms_int)



cdef class Box:
    """Represents the shape of the simulation cell

    Parameters
    ----------
    lammps : :class:`lammps.core.Lammps`
         LAMMPS Object
    """
    cdef Lammps lammps

    def __cinit__(self, Lammps lammps):
        """ Docstring in Box base class (sphinx can find doc when compiled) """
        self.lammps = lammps

    def from_lattice_const(self, atom_types, lengths, angles=None, region_id='lammps-box'):
        """ Create Unit cell (can only be run once)

        Parameters
        ----------
        atom_types : int
            number of atom types for simulation (I know weird place to require)
        lengths : vector
            lengths length 3 (a, b, c) lattice constants
        angles : vector
            angles length 3 (alpha, beta, gamma) lattice constants default (pi/2, pi/2, pi/2)
        region_id : str
            name of region for lammps simulation. default 'lammps-box'

        Essentially runs:
        region region-id prism xlo xhi ylo yhi zlo zhi xy xz yz
        create_box N region-ID keyword value ...

..      todo::
        This functions needs more thought. It does too much but within
        the LAMMPS C++ function calls I can't see a nice way to do it.
        """
        if angles is None:
            angles = [pi/2., pi/2., pi/2.]

        (lx, ly, lz), (xy, xz, yz) = lattice_const_to_lammps_box(lengths, angles)
        self.lammps.command((
            "region {} prism 0.0 {} 0.0 {} 0.0 {} {} {} {}"
        ).format(region_id, lx, ly, lz, xy, xz, yz))

        self.lammps.command((
            "create_box {} {}"
        ).format(atom_types, region_id))

    def update_lattice_const(self, lengths, angles=None):
        """Update Unit cell

        Parameters
        ----------
        lengths : vector
            lengths length 3 (a, b, c) lattice constants
        angles : vector
            angles length 3 (alpha, beta, gamma) lattice constants default (pi/2, pi/2, pi/2)

        Does not adjust the atom coordinates. However it does set the
        global box, set processor grid, and set local box correctly.

.. warning::

        This does not update the atom coordinates. This means that they can easily be lost!
        """
        if angles == None:
            angles = [pi/2., pi/2., pi/2.]

        (lx, ly, lz), (xy, xz, yz) = lattice_const_to_lammps_box(lengths, angles)
        cdef double* boxlo = [0., 0., 0.]
        cdef double* boxhi = [lx, ly, lz]
        lammps_reset_box(self.lammps._lammps, boxlo, boxhi, xy, xz, yz)

    @property
    def dimension(self):
        """ The dimension of the lammps run either 2D or 3D"""
        return self.lammps._lammps.domain.dimension

    @property
    def lohi(self):
        """ LAMMPS box description of boxhi and boxlo

        Additional documentation can be found at `lammps
        <http://lammps.sandia.gov/doc/Section_howto.html?highlight=prism#howto-12>`_. For
        example return numpy array would be

..      code-block:: python

        lohi = np.array([
            [0.0, 0.0, 0.0],
            [10.0, 10.0, 10.0]
        ])
        """
        cdef int dim = self.dimension
        cdef double[::1] boxlo = <double[:dim]>self.lammps._lammps.domain.boxlo
        cdef double[::1] boxhi = <double[:dim]>self.lammps._lammps.domain.boxhi
        return np.array([boxlo, boxhi])

    @property
    def tilts(self):
        """ LAMMPS box description of xy xz yz

        Additional documentation can be found at `lammps
        <http://lammps.sandia.gov/doc/Section_howto.html?highlight=prism#howto-12>`_.
        For example a return dictionary would be

..      code-block:: python

        tilts = np.array([0.0, 0.0, 0.0])

..      todo::

        how to handle 2d? Needs to be addressed throughout
        """
        cdef double xy = self.lammps._lammps.domain.xy
        cdef double xz = self.lammps._lammps.domain.xz
        cdef double yz = self.lammps._lammps.domain.yz
        return np.array([xy, xz, yz])

    @property
    def lengths_angles(self):
        """ Calculates the lengths and angles from lammps box

        Additional documentation can be found at `lammps
        <http://lammps.sandia.gov/doc/Section_howto.html?highlight=prism#howto-12>`_. See
        implementation :py:func:`lammps_box_to_lattice_const`.
        """
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

    @property
    def lengths(self):
        """Calculates lengths from lammps box. See :py:func:`lengths_angles`"""
        return self.lengths_angles[0]

    @property
    def angles(self):
        """ Calculates lengths from lammps box. See :py:func:`lengths_angles`"""
        return self.lengths_angles[1]

    @property
    def volume(self):
        """ Volume of box given in lammps units"""
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
    cdef char **argv = <char**>malloc(len(args) * sizeof(char*))
    cdef char* temp
    cdef int i
    for i in range(len(args)):
        byte_string = args[i].encode('UTF-8')
        temp = byte_string
        argv[i] = <char*>malloc(strlen(temp) + 1)
        strcpy(argv[i], temp)
    return argv
