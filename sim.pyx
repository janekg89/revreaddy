# distutils: language = c++
# distutils: extra_compile_args = -std=c++11
# distutils: sources = src/Observable.cpp src/Simulation.cpp src/Particle.cpp src/Random.cpp src/Force.cpp src/Trajectory.cpp src/RadialDistribution.cpp src/utils.cpp
# distutils: include_dirs = /usr/include /usr/local/include include
# distutils: libraries = m gsl gslcblas

from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp cimport bool
cimport numpy as np

cdef extern from "Simulation.h":
	cdef cppclass Simulation:
		Simulation() except +
		unsigned long int maxTime
		double temperature
		double timestep
		bool isPeriodic
		double boxsize
		double repulsionStrength
		void addParticle(vector[double], string, double, double)
		void run()
		vector[double] getPosition(int)
		int getParticleNumber()
		void pushObservable(Observable)

cdef class pySimulation:
	cdef Simulation *thisptr
	def __cinit__(self):
		self.thisptr = new Simulation()
	def __dealloc__(self):
		del self.thisptr
	def addParticle(
		self,
		initialPosition,
		particleType,
		radius,
		diffusionConst):
		self.thisptr.addParticle(
			initialPosition,
			particleType,
			radius,
			diffusionConst)
	def run(self):
		self.thisptr.run()
	def getPosition(self, index):
		return self.thisptr.getPosition(index)
	def getParticleNumber(self):
		return self.thisptr.getParticleNumber()
	#def pushObservable(self, Observable observable):
	#	self.thisptr.pushObservable(observable)
	property maxTime:
		def __get__(self): return self.thisptr.maxTime
		def __set__(self, maxTime): self.thisptr.maxTime = maxTime
	property temperature:
		def __get__(self): return self.thisptr.temperature
		def __set__(self, temperature): self.thisptr.temperature = temperature
	property timestep:
		def __get__(self): return self.thisptr.timestep
		def __set__(self, timestep): self.thisptr.timestep = timestep
	property isPeriodic:
		def __get__(self): return self.thisptr.isPeriodic
		def __set__(self, isPeriodic): self.thisptr.isPeriodic = isPeriodic
	property boxsize:
		def __get__(self): return self.thisptr.boxsize
		def __set__(self, boxsize): self.thisptr.boxsize = boxsize
	property repulsionStrength:
		def __get__(self): return self.thisptr.repulsionStrength
		def __set__(self, repulsionStrength): 
			self.thisptr.repulsionStrength = repulsionStrength

#the following does not yet work
cdef extern from "Observable.h":
	cdef cppclass Observable:
		Observable() except +

cdef extern from "Trajectory.h":
	cdef cppclass Trajectory(Observable):
		Trajectory() except +
		void writeBufferToFile()

cdef class pyObservable:
	pass

cdef class pyTrajectory:
	cdef Trajectory *thisptr
	def __cinit__(self):
		self.thisptr = new Trajectory()
	def __dealloc__(self):
		del self.thisptr
	def writeBufferToFile(self):
		self.thisptr.writeBufferToFile()

