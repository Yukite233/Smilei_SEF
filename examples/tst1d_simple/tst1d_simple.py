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
Lsim = [200.*l0]	# length of the simulation
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
dim = '1d3v'

# order of interpolation
#
interpolation_order = 1

# SIMULATION BOX : for all space directions (use vector)
# cell_length: length of the cell
# sim_length: length of the simulation in units of the normalization wavelength 
#
cell_length = [1]
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
bc_em_type_x = ['Dirichlet', 'Dirichlet']
bc_em_value_x = [0.0, 0.0]


#Topology:
#number_of_procs: Number of MPI processes in each direction.
#clrw: width of a cluster in number of cell. Warning: clrw must divide nspace_win_x.
number_of_procs = [4]


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

