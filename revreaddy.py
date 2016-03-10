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

	# config properties
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
	def kt(self, kt):
		if (kt <= 0.):
			raise Exception("The temperature must be positive.")
		self.config.setKT(kt)

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
	def boxsize(self, bs):
		if (bs <= 0.):
			raise Exception("The boxsize must be positive.")
		self.config.setBoxsize(bs)
	
	# simulation properties