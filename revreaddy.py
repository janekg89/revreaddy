"""
A particle-based reaction-diffusion simulation.

This module is the convenience/user layer of the underlying
C++ revreaddy library (C-extension). It offers a simple 
interface and useful methods to run reaction-diffusion
simulations, save their results in compact binary files
and save or load earlier states of a simulation with ease.

It defines the main class Sim, which communicates with
the C-extension, and methods that act on Sim that are
considered 'utilities'.

Dependencies are:
* revreaddyPy (C-extension)
* numpy
* h5py

"""
import numpy as np
import h5py
import revreaddyPy as ext
import logging
# this is only the log-level of the python module
# the logging level for the C-extension is set on compile time
logging.basicConfig(
	format='%(asctime)s %(levelname)s: %(message)s',
	datefmt='%Y/%m/%d %I:%M:%S',
	level=logging.INFO)
import time

class Sim(object):
	"""
	User interface class to setup and run simulations.

	Main object that wraps the three most important objects
	of the C-extension: Config, World and Simulation.
	Config represents the data that does not change during
	the course of the simulation (particle types, 
	reactions, interactions, temperature ...). World is the
	data that does change during the simulation like
	the position of particles and the cumulative runtime.
	Both Config and World represent a particular system.
	Simulation contains the methods to execute the 
	simulation. There are different integrators that can
	be chosen.

	"""
	def __init__(self, which_impl=""):
		"""Construct the three main objects and choose implementation"""
		self.config = ext.Config()
		self.world = ext.World()
		self.simulation = ext.Simulation(self.world, self.config, which_impl)

	# wrapped config properties
	@property
	def timestep(self):
	    return self.config.getTimestep()
	@timestep.setter
	def timestep(self, dt):
		if (dt <= 0.):
			raise Exception("Timestep must be positive.")
		self.config.setTimestep(dt)

	@property
	def kt(self):
	    return self.config.getKT()
	@kt.setter
	def kt(self, kt_):
		if (kt_ <= 0.):
			raise Exception("The temperature must be positive.")
		self.config.setKT(kt_)

	@property
	def is_periodic(self):
	    return self.config.getIsPeriodic()
	@is_periodic.setter
	def is_periodic(self, periodic):
		self.config.setIsPeriodic(periodic)

	@property
	def boxsize(self):
	    return self.config.getBoxsize()
	@boxsize.setter
	def boxsize(self, boxsize_):
		if (boxsize_ <= 0.):
			raise Exception("The boxsize must be positive.")
		self.config.setBoxsize(boxsize_)
	
	# wrapped simulation properties
	@property
	def use_neighborlist(self):
	    return self.simulation.getUseNeighborlist()
	@use_neighborlist.setter
	def use_neighborlist(self, nl):
		self.simulation.setUseNeighborlist(nl)

	# wrapped config methods
	def delete_all_particle_types(self):
		self.config.deleteAllParticleTypes()

	def new_type(self, name, radius, diffusion_constant):
		self.config.new_Type(name, radius, diffusion_constant)

	def delete_all_geometries(self):
		self.config.deleteAllGeometries()

	def new_wall(self, normal, point, strength, particle_type_ids):
		normal = np.array(normal)
		if (normal.shape != (3,)):
			raise Exception("Normal vector has wrong shape.")
		point = np.array(point)
		if (point.shape != (3,)):
			raise Exception("Point vector has wrong shape.")
		# use python built-in type int because boost doesn't know numpy types
		particle_type_ids = np.array(particle_type_ids, dtype=int)
		if ( len(particle_type_ids.shape) != 1):
			raise Exception("Particle types must be a one-dimensional container.")
		self.config.new_Wall(normal, point, strength, particle_type_ids)

	def new_double_well_z(self, distance_minima, strength, particle_type_ids):
		particle_type_ids = np.array(particle_type_ids, dtype=int)
		if ( len(particle_type_ids.shape) != 1):
			raise Exception("Particle types must be a one-dimensional container.")
		self.config.new_DoubleWellZ(distance_minima, strength, particle_type_ids)

	def delete_all_interactions(self):
		self.config.deleteAllInteractions()

	def new_soft_repulsion(self, name, affected_tuple, repulsion_strength):
		affected_tuple = np.array(affected_tuple, dtype=int)
		if (affected_tuple.shape != (2,)):
			raise Exception("Affected tuple has wrong shape.")
		self.config.new_SoftRepulsion(name, affected_tuple, repulsion_strength)

	def new_lennard_jones(self, name, affected_tuple, epsilon):
		affected_tuple = np.array(affected_tuple, dtype=int)
		if (affected_tuple.shape != (2,)):
			raise Exception("Affected tuple has wrong shape.")
		self.config.new_LennardJones(name, affected_tuple, epsilon)

	def delete_all_reactions(self):
		self.config.deleteAllReactions()

	def new_conversion(self, name, forward_type, backward_type, forward_rate, backward_rate):
		"""Register a new conversion reaction to the config."""
		self.config.new_Conversion(name, forward_type, backward_type, forward_rate, backward_rate)

	def new_fusion(self, name, forward_type_a, forward_type_b, backward_type_c,	forward_rate, backward_rate, reaction_distance):
		"""Register a new fusion reaction to the config."""
		self.config.new_fusion(name, forward_type_a, forward_type_b, backward_type_c, forward_rate, backward_rate, reaction_distance)

	# TODO WIP this only exists as long as fusion is configured manually
	def configure_fusion(self, reaction_index, interaction_indices, inverse_partition, max_distr, radii_sum, reaction_radii_sum, mean_distr, inverse_temperature, radius_a, radius_b):
		interaction_indices = np.array(interaction_indices, dtype=int)
		if (len(interaction_indices.shape) != 1):
			raise Exception("Interaction-indices must be a one-dimensional container.")
		self.config.configureFusion(reaction_index, interaction_indices, inverse_partition, max_distr, radii_sum, reaction_radii_sum, mean_distr, inverse_temperature, radius_a, radius_b)

	# wrapped world methods
	def delete_all_particles(self):
		self.world.deleteAllParticles()

	def add_particle(self,init_pos, particle_type_id):
		init_pos = np.array(init_pos)
		if (init_pos.shape != (3,)):
			raise Exception("Initial position has the wrong shape.")
		self.world.addParticle(init_pos, particle_type_id)

	# wrapped simulation methods
	def run(self, number_timesteps, timestep_=None):
		"""Start the simulation."""
		if (timestep_ is not None):
			self.timestep = timestep_
		if (number_timesteps <= 0):
			raise Exception("Number of timesteps must be positive.")
		logging.info("Run with timestep "+str(self.timestep)+" and "+str(number_timesteps)+" timesteps")
		t1 = time.clock()
		self.simulation.run(number_timesteps)
		t2 = time.clock()
		logging.info("Finished after "+str(t2-t1)+" seconds.")

	def delete_all_observables(self):
		self.simulation.deleteAllObservables()

	# TODO WIP this will become obsolete when write-periods are implemented
	def write_observables_to_file(self):
		self.simulation.writeAllObservablesToFile()

	def new_trajectory(self, rec_period, filename):
		self.simulation.new_Trajectory(rec_period, filename)

	def new_radial_distribution(self, rec_period, filename, ranges, considered):
		"""
		Register an observable that calculates the radial distribution function.

		'considered' contains a list of pairs of particle types, that are
		considered in the calculation of the RDF. Here the order or types
		is significant. I.e. if considered is [[0,1]], then in the loop over
		all particles only (0,1) pairs are considered but not (1,0) pairs.
		Hence when 'considered' is [[0,0]] the distances between a 
		particular pair of particles with types==0 is counted twice.

		"""
		ranges = np.array(ranges)
		if (len(ranges.shape) != 1):
			raise Exception("Ranges must be a one-dimensional container.")	
		considered = np.array(considered, dtype=int)
		if ( (len(considered.shape) != 2) | (len(considered.shape[1] != 2)) ):
			raise Exception("Considered types must have a shape of (n,2).")
		self.simulation.new_RadialDistribution(rec_period, filename, ranges, considered)

	def new_mean_squared_displacement(self, rec_period, filename, particle_type_id):
		"""
		Register an observable that calculates the mean squared displacement.

		After the simulation has started this observable keeps track
		of all particles that are of type 'particle_type_id' initially.
		If they are deleted they will be ignored further on, effectively
		decreasing the number of considered particles over time. If
		new particles of type 'particle_type_id' are created during
		runtime they are ignored.
		When using periodic boundary conditions, the 'real' traveled
		distance is calculated, so you don't need to worry about
		boundary effects here. 

		"""
		self.simulation.new_MeanSquaredDisplacement(rec_period, filename, particle_type_id)

	def new_probability_density(self, rec_period, filename, particle_type_id, ranges, coord):
		"""
		Register an observable that records a histogram of particle-positions.

		This simply computes a histogram of one coordinate 'coord' (0, 1 or 2)
		of particle-positions that have the particle type 'particle_type_id'.
		At every recording event, all particles of this type are considered even
		when the number of particles changes. 

		"""
		ranges = np.array(ranges)
		if (len(ranges.shape) != 1):
			raise Exception("Ranges must be a one-dimensional container.") 
		if (np.sort(ranges) == ranges):
			raise Exception("Ranges must be sorted/monotonically increasing.")
		self.simulation.new_ProbabilityDensity(rec_period, filename, particle_type_id, ranges, coord)

	def new_energy(self, rec_period, filename):
		"""Register an observable that records the total internal energy."""
		self.simulation.new_Energy(rec_period, filename)

	def new_acceptance(self, rec_period, filename, reactions_or_diffusion):
		"""
		Register an observable that records the acceptance probability.

		Note that this keeps track of the calculated acceptance value
		not the real acceptance that would have to be averaged over 
		many timesteps. 'reactions_or_diffusion' decides if
		the acceptance of the reaction or the diffusion step is
		considered.

		"""
		self.simulation.new_Acceptance(rec_period, filename, reactions_or_diffusion)

	def new_particle_numbers(self, rec_period, filename, particle_type_id):
		"""Register an observable that records the number of particles."""
		self.simulation.new_ParticleNumbers(rec_period, filename, particle_type_id)

	# derived methods
	def show_config(self):
		pass