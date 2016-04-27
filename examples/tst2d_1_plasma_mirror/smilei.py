"""@package pyinit
    Definition of Smilei components
"""


class SmileiComponentType(type):
    """Metaclass to all Smilei components"""
    
    # Constructor of classes
    def __init__(self, name, bases, attrs):
        self._list = []
        self.verify = True
        self.current = 0
        # Run standard metaclass init
        super(SmileiComponentType, self).__init__(name, bases, attrs)
    
    # Functions to define the iterator
    def __iter__(self):
        return self
    def next(self):
        if self.current >= len(self._list):
            raise StopIteration
        self.current += 1
        return self._list[self.current - 1]
    
    # Function to return one given instance, for example DiagParticles[0]
    # Special case: species can also be indexed by their name: Species["ion1"]
    def __getitem__(self, key):
        if self.__name__ == "Species" and type(key) is str:
            for obj in self._list:
                if obj.species_type == key:
                    return obj
        else:
            return self._list[key]
    
    # Function to return the number of instances, for example len(Species)
    def __len__(self):
        return len(self._list)
    
    # Function to display the content of the component
    def __repr__(self):
        if len(self._list)==0:
            return "<Empty list of "+self.__name__+">"
        else:
            l = []
            for obj in self._list: l.append(str(obj))
            return "["+", ".join(l)+"]"


class SmileiComponent(object):
    """Smilei component generic class"""
    __metaclass__ = SmileiComponentType
    
    # This constructor is used always for all child classes
    def __init__(self, **kwargs):
        if kwargs is not None: # add all kwargs as internal class variables
            for key, value in kwargs.iteritems():
                if key=="_list":
                    print "Python warning: in "+type(self).__name__+": cannot have argument named '_list'. Discarding."
                else:
                    setattr(self, key, value)
        type(self)._list.append(self) # add the current object to the static list "list"


class Species(SmileiComponent):
    """Species parameters"""
    species_type = None
    initPosition_type = None
    initMomentum_type = ""
    n_part_per_cell = None
    c_part_max = 1.0
    charge_density = None
    nb_density = None
    mean_velocity = [0.]
    temperature = [1e-10]
    dynamics_type = "norm"
    time_frozen = 0.0
    radiating = False
    bc_part_type_west = None
    bc_part_type_east = None
    bc_part_type_north = None
    bc_part_type_south = None
    ionization_model = "none"
    atomic_number = None
    isTest = False


class Laser(SmileiComponent):
    """Laser parameters"""
    boxSide = None
    a0 = None
    omega0 = 1.
    delta = 1.
    tchirp = 0.
    focus = []
    angle = 0.
    delay = 0.
    time_profile = None
    int_params = []
    double_params = []
    transv_profile = None
    int_params_transv = []
    double_params_transv = []

class Collisions(SmileiComponent):
    """Collisions parameters"""
    species1 = None
    species2 = None
    coulomb_log = 0.
    debug_every = 0


#diagnostics
class DiagProbe(SmileiComponent):
    """Diagnostic probe"""
    every = 0
    time_range = [None, None]
    number = []
    pos = []
    pos_first = []
    pos_second = []
    pos_third = []
    fields = []

class DiagParticles(SmileiComponent):
    """Diagnostic particles"""
    output = None
    every = None
    time_average = 1
    species = None
    axes = []

class DiagPhase(SmileiComponent):
    """Diagnostic phase"""
    pass

class DiagScalar(SmileiComponent):
    """Diagnostic scalar"""
    every = None
    time_range = []
    precision = 10
    vars = []

# external fields
class ExtField(SmileiComponent):
    """External Field"""
    field = []
    profile = None

# default simulation values
output_script = "smilei.py"
dump_step = 0
dump_minutes = 0.0
exit_after_dump = True
restart = False
check_stop_file = False
dump_file_sequence = 2
sim_units = ""
wavelength_SI = 0.
dim = ""
interpolation_order = 2
res_time = None
res_space = []
timestep = None
cell_length = []
sim_time = None
sim_length = []
bc_em_type_x = []
bc_em_type_y = []
time_fields_frozen = 0.0
nspace_win_x = 0
t_move_win = 0.0
vx_win = 1.
clrw = 1
every = 0
number_of_procs = [None]
print_every = None
fieldDump_every = 0
fieldsToDump = []
avgfieldDump_every = None
ntime_step_avg = 0
particleDump_every = None # for backwards-compatibility
time_fields_frozen = 0.
# Some predefined profiles (see doc)

def constant(value, xvacuum=0., yvacuum=0.):
    global dim, sim_length
    if dim == "1d3v": return lambda x: value if x>=xvacuum else 0.
    if dim == "2d3v": return lambda x,y: value if (x>=xvacuum and y>=yvacuum) else 0.

def trapezoidal(max,
                xvacuum=0., xplateau=None, xslope1=0., xslope2=0.,
                yvacuum=0., yplateau=None, yslope1=0., yslope2=0. ):
    global dim, sim_length
    if len(sim_length)>0 and xplateau is None: xplateau = sim_length[0]-xvacuum
    if len(sim_length)>1 and yplateau is None: yplateau = sim_length[1]-yvacuum
    def fx(x):
        # vacuum region
        if x < xvacuum: return 0.
        # linearly increasing density
        elif x < xvacuum+xslope1: return max*(x-xvacuum) / xslope1
        # density plateau
        elif x < xvacuum+xslope1+xplateau: return max
        # linearly decreasing density
        elif x < xvacuum+xslope1+xplateau+xslope2:
            return max*(1. - ( x - (xvacuum+xslope1+xslope2) ) / xslope2)
        # beyond the plasma
        else: return 0.0
    def fy(y):
        # vacuum region
        if y < yvacuum: return 0.
        # linearly increasing density
        elif y < yvacuum+yslope1: return (y-yvacuum) / yslope1
        # density plateau
        elif y < yvacuum+yslope1+yplateau: return 1.
        # linearly decreasing density
        elif y < yvacuum+yslope1+yplateau+yslope2:
            return 1. - ( y - (yvacuum+yslope1+yslope2) ) / yslope2
        # beyond
        else: return 0.0
    if dim == "1d3v": return fx
    if dim == "2d3v": return lambda x,y: fx(x)*fy(y)

def gaussian(max,
             xvacuum=0., xlength=None, xfwhm=None, xcenter=None, xorder=2,
             yvacuum=0., ylength=None, yfwhm=None, ycenter=None, yorder=2 ):
    import math
    global dim, sim_length
    if len(sim_length)>0:
        if xlength is None: xlength = sim_length[0]-xvacuum
        if xfwhm   is None: xfwhm   = xlength/3.
        if xcenter is None: xcenter = xvacuum + xlength/2.
    if len(sim_length)>1: 
        if ylength is None:ylength = sim_length[1]-yvacuum
        if yfwhm   is None:yfwhm   = ylength/3.
        if ycenter is None:ycenter = yvacuum + ylength/2.
    def fx(x):
        sigmaN = (0.5*xfwhm)**xorder/math.log(2.0)
        # vacuum region
        if x < xvacuum: return 0.
        # gaussian
        elif x < xvacuum+xlength: return max*math.exp( -(x-xcenter)**xorder / sigmaN )
        # beyond
        else: return 0.0
    def fy(y):
        if yorder == 0: return 1.
        sigmaN = (0.5*yfwhm)**yorder/math.log(2.0)
        # vacuum region
        if y < yvacuum: return 0.
        # gaussian
        elif y < yvacuum+ylength: return math.exp( -(y-ycenter)**yorder / sigmaN )
        # beyond
        else: return 0.0
    if dim == "1d3v": return fx
    if dim == "2d3v": return lambda x,y: fx(x)*fy(y)

def polygonal(xpoints=[], xvalues=[]):
    global dim, sim_length
    if len(xpoints)!=len(xvalues):
        raise Exception("polygonal profile requires as many points as values")
    if len(sim_length)>0 and len(xpoints)==0:
        xpoints = [0., sim_length[0]]
        xvalues = [1., 1.]
    N = len(xpoints)
    def f(x,y=0.):
        # vacuum region
        if x < xpoints[0]: return 0.0;
        # polygon region (defined over N segments)
        elif x < xpoints[-1]:
            for i in range(1,len(xpoints)):
                if xpoints[i-1]==xpoints[i]: return xvalues[i-1]
                if x>=xpoints[i-1] and x<xpoints[i]:
                    m = (xvalues[i]-xvalues[i-1])/(xpoints[i]-xpoints[i-1])
                    return xvalues[i-1] + m * ( x-xpoints[i-1] )
        # beyond
        else: return 0.
    return f

def cosine(base, amplitude=1.,
           xvacuum=0., xlength=None, phi=0., xnumber=1):
    import math
    global sim_length
    if len(sim_length)>0 and xlength is None: xlength = sim_length[0]-xvacuum
    def f(x,y=0.):
        #vacuum region
        if x < xvacuum: return 0.
        # profile region
        elif x < xvacuum+xlength:
            return base + amplitude * math.cos(phi + 2.*math.pi * xnumber * (x-xvacuum)/xlength)
        # beyond
        else: return 0.
    return f
############### BEGIN USER NAMELISTS ###############
# ----------------------------------------------------------------------------------------
# 					SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------
#
# Remember: never override the following names:
#           SmileiComponent, Species, Laser, Collisions, DiagProbe, DiagParticles,
#           DiagScalar, DiagPhase or ExtField
#
import math

l0 = 1.0			# normalized by 1.0e-5 m
t0 = 1.0e-12 * 3.0e8 / 1.0e-5	# real time = 1.0e-12 s
Lsim = [50.*l0,20.*l0]	# length of the simulation
Tsim = 1000.*t0			# duration of the simulation
resx = 20.				# nb of cells in on laser wavelength
rest = 30.				# time of timestep in one optical cycle 

wavelength_SI = 1.e-6


# dim: Geometry of the simulation
#      1d3v = cartesian grid with 1d in space + 3d in velocity
#      2d3v = cartesian grid with 2d in space + 3d in velocity
#      3d3v = cartesian grid with 3d in space + 3d in velocity
#      2drz = cylindrical (r,z) grid with 3d3v particles
#
dim = '2d3v'

# order of interpolation
#
interpolation_order = 1

# SIMULATION BOX : for all space directions (use vector)
# cell_length: length of the cell
# sim_length: length of the simulation in units of the normalization wavelength 
#
cell_length = [1,1]
sim_length  = Lsim

# SIMULATION TIME
# timestep: duration of the timestep
# sim_time: duration of the simulation in units of the normalization period 
#
timestep = t0
sim_time = Tsim
 
# ELECTROMAGNETIC BOUNDARY CONDITIONS
# bc_em_type_x/y/z : boundary conditions used for EM fields 
#                    periodic = periodic BC (using MPI topology)
#                    silver-muller = injecting/absorbing BC
#                    reflective = consider the ghost-cells as a perfect conductor
#
bc_em_type_x = ['silver-muller']
bc_em_type_y = ['silver-muller']


#Topology:
#number_of_procs: Number of MPI processes in each direction.
#clrw: width of a cluster in number of cell. Warning: clrw must divide nspace_win_x.
number_of_procs = [2, 2]


# RANDOM seed 
# this is used to randomize the random number generator
random_seed = 0



# DEFINE ALL SPECIES
# species_type       = string, given name to the species (e.g. ion, electron, positron, test ...)
# initPosition_type  = string, "regular" or "random"
# initMomentum_type  = string "cold", "maxwell-juettner" or "rectangular"
# c_part_max         = float, factor on the memory reserved for the total number of particles
# mass               = float, particle mass in units of the electron mass
# dynamics_type      = string, type of species dynamics = "norm" or "rrLL"
# time_frozen        = float, time during which particles are frozen in units of the normalization time
# radiating          = boolean, if true, incoherent radiation calculated using the Larmor formula 
# n_part_per_cell    = integer or function, number of particles/cell
# charge             = float or function, particle charge in units of the electron charge
# charge_density     = float or function, species charge density in units of the "critical" density
#     or nb_density for number density
# mean_velocity      = list of floats or functions, mean velocity in units of the speed of light
# temperature        = list of floats or functions, temperature in units of m_e c^2
# Predefined functions: constant, trapezoidal, gaussian, polygonal, cosine
#
Species(
	species_type = 'ion',
	initPosition_type = 'random',
	initMomentum_type = 'maxwell-juettner',
	ionization_model = 'none',
	n_part_per_cell = 80,
	c_part_max = 1.0,
	mass = 1836.0,
	charge = 1.0,
	nb_density = 3.5e-5,
	temperature = [1.9e-5],
	time_frozen = 0.0,
	bc_part_type_west  = 'supp',
	bc_part_type_east  = 'supp',
	bc_part_type_south = 'refl',
	bc_part_type_north = 'refl'
)
Species(
	species_type = 'eon',
	initPosition_type = 'random',
	initMomentum_type = 'maxwell-juettner',
	ionization_model = 'none',
	n_part_per_cell = 80,
	c_part_max = 1.0,
	mass = 1.0,
	charge = -1.0,
	nb_density = 3.5e-5,
	temperature = [1.9e-5],
	time_frozen = 0.,
	bc_part_type_west  = 'supp',
	bc_part_type_east  = 'supp',
	bc_part_type_south = 'refl',
	bc_part_type_north = 'refl'
)


Species(
	species_type = 'neutral',
	initPosition_type = 'random',
	initMomentum_type = 'maxwell-juettner',
	ionization_model = 'none',
	n_part_per_cell = 80,
	c_part_max = 1.0,
	mass = 1836.0,
	charge = 0.0,
	nb_density = 3.5e-5,
	temperature = [1.9e-5],
	time_frozen = 0.0,
	bc_part_type_west  = 'supp',
	bc_part_type_east  = 'supp',
	bc_part_type_south = 'refl',
	bc_part_type_north = 'refl'
)




# COLLISIONS
# species1    = list of strings, the names of the first species that collide
# species2    = list of strings, the names of the second species that collide
#               (can be the same as species1)
# coulomb_log = float, Coulomb logarithm. If negative or zero, then automatically computed.
Collisions(
	species1 = ["eon"],
	species2 = ["ion"],
	coulomb_log = 5,
	collisions_type = "coulomb"
)
Collisions(
	species1 = ["eon"],
	species2 = ["eon"],
	coulomb_log = 1,
	collisions_type = "coulomb"
)
Collisions(
	species1 = ["ion"],
	species2 = ["ion"],
	coulomb_log = 1,
	collisions_type = "coulomb"
)



# ---------------------
# DIAGNOSTIC PARAMETERS
# ---------------------

# print_every (on screen text output) 
# print_every = 60

dump_step = 100


################ END USER NAMELISTS ################
"""@package pycontrol
    here we check if the namelist is clean and without errors
"""

import gc 
gc.collect()

def _smilei_check():
    """Do checks over the script"""
    
    # Verify classes were not overriden
    for CheckClassName,CheckClass in {"SmileiComponent":SmileiComponent,"Species":Species,
            "Laser":Laser,"Collisions":Collisions,"DiagProbe":DiagProbe,"DiagParticles":DiagParticles,
            "DiagScalar":DiagScalar,"DiagPhase":DiagPhase,"ExtField":ExtField}.iteritems():
        try:
            if not CheckClass.verify: raise
        except:
            raise Exception("ERROR in the namelist: it seems that the name `"+CheckClassName+"` has been overriden")
    
    # Check species for undefined/duplicate species_type
    all_species=[]
    for spec in Species:
        if spec.species_type == None:
            raise Exception("ERROR in the namelist: there is a species without species_type")
        elif spec.species_type in all_species:
            raise Exception("ERROR in the namelist: there is duplicate species_type")
        else:
            all_species.append(spec.species_type)
    
# this function will be called after initialising the simulation, just before entering the time loop
# if it returns false, the code will call a Py_Finalize();
def _keep_python_running():
    for las in Laser:
        for prof in (las.time_profile, las.transv_profile):
            if callable(prof): return True
    return False

# Prevent creating new components (by mistake)
def _noNewComponents(cls, *args, **kwargs):
    print "Please do not create a new "+cls.__name__
    return None
SmileiComponent.__new__ = staticmethod(_noNewComponents)


