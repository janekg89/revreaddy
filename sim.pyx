# distutils: language = c++
# distutils: extra_compile_args = -std=c++11
# distutils: sources = src/Observable.cpp src/Simulation.cpp src/Particle.cpp src/Random.cpp src/Force.cpp src/Trajectory.cpp src/TrajectorySingle.cpp src/RadialDistribution.cpp src/utils.cpp
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
		bool verbose
		unsigned long int acceptions
		unsigned long int rejections

		void addParticle(vector[double], string, double, double)
		void run()
		vector[double] getPosition(int)
		void setPosition(int, vector[double])
		int getParticleNumber()
		void deleteAllParticles()
		void writeAllObservablesToFile()
		string showObservables()
		void deleteAllObservables()
		void new_Trajectory(string)
		void new_TrajectorySingle()
		vector[vector[double]] getTrajectorySingle()

cdef class pySimulation:
	cdef Simulation *thisptr
	def __cinit__(self):
		self.thisptr = new Simulation()
	def __dealloc__(self):
		del self.thisptr
	def addParticle(
		self,
		initialPosition,
		particleType='soft',
		radius=1.,
		diffusionConst=1.):
		self.thisptr.addParticle(
			initialPosition,
			particleType,
			radius,
			diffusionConst)
	def run(self):
		self.thisptr.run()
	def getPosition(self, index): 
		return self.thisptr.getPosition(index)
	def setPosition(self, index, newPos):
		self.thisptr.setPosition(index, newPos)
	def getParticleNumber(self): 
		return self.thisptr.getParticleNumber()
	def deleteAllParticles(self):
		self.thisptr.deleteAllParticles()
	def writeAllObservablesToFile(self): 
		self.thisptr.writeAllObservablesToFile()
	def showObservables(self):
		return self.thisptr.showObservables()
	def deleteAllObservables(self): 
		self.thisptr.deleteAllObservables()
	def new_Trajectory(self, filename): 
		self.thisptr.new_Trajectory(filename)
	def new_TrajectorySingle(self):
		self.thisptr.new_TrajectorySingle()
	def getTrajectorySingle(self):
		return self.thisptr.getTrajectorySingle()
		
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
	property verbose:
		def __get__(self): return self.thisptr.verbose
		def __set__(self, verbose): self.thisptr.verbose = verbose
	property acceptions:
		def __get__(self): return self.thisptr.acceptions
		def __set__(self, acceptions): self.thisptr.acceptions = acceptions
	property rejections:
		def __get__(self): return self.thisptr.rejections
		def __set__(self, rejections): self.thisptr.rejections = rejections

	# derived functions
	def acceptanceRate(self):
		acc = 1./(1.+ float(self.rejections)/float(self.acceptions) )
		return round(acc, 5)
