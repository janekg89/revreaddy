#!python
#cython: embedsignature=True
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp cimport bool
cimport numpy as np
import time

cdef extern from "Simulation.h":

	cdef cppclass World:
		World() except +

		unsigned long int acceptionsDynamics
		unsigned long int rejectionsDynamics
		unsigned long int acceptionsReactions
		unsigned long int rejectionsReactions
		double cumulativeRuntime

		void addParticle(vector[double], unsigned long int)
		void removeParticle(unsigned long int)
		vector[double] getPosition(unsigned long int)
		void setPosition(unsigned long int, vector[double])
		unsigned int getTypeId(unsigned long int)
		void setTypeId(unsigned long int, unsigned int)
		void deleteAllParticles()

	cdef cppclass Config:
		Config(World*) except +

		unsigned long int maxTime
		double temperature
		double kBoltzmann
		double timestep
		bool isPeriodic
		double boxsize
		bool isReversibleDynamics
		bool isReversibleReactions
		bool useNeighborList
		unsigned int reactionPropagation

		void new_Type(string, double, double, double)
		unsigned int getNumberOfTypes()
		string getDictName(unsigned int)
		double getDictRadius(unsigned int)
		double getDictDiffusionConstant(unsigned int)
		double getDictReactionRadius(unsigned int)
		unsigned int getParticleNumber()
		void writeAllObservablesToFile()
		void writeLastObservableToFile()
		string showObservables()
		void deleteAllObservables()
		void deleteLastObservable()
		void new_Trajectory(unsigned long int, string)
		void new_RadialDistribution(
			unsigned long int,
			string,
			vector[double],
			vector[vector[uint]])
		void new_MeanSquaredDisplacement(
			unsigned long int, string, unsigned int)
		void new_ProbabilityDensity(
			unsigned long int,
			string,
			unsigned int,
			vector[double],
			unsigned int)
		void new_Energy(unsigned long int, string)
		void new_Acceptance(unsigned long int, string)
		void deleteAllGeometries()
		void new_Wall(
			vector[double], vector[double],
			double, vector[uint])
		void new_DoubleWellZ(double, double, vector[uint])
		void deleteAllForces()
		void new_SoftRepulsion(string, vector[uint], double)
		void new_LennardJones(string, vector[uint], double)
		unsigned int getNumberForces()
		string getForceName(unsigned int)
		string getForceType(unsigned int)
		vector[uint] getForceAffectedTuple(unsigned int)
		vector[double] getForceParameters(unsigned int)

	cdef cppclass Simulation:
		Simulation(World*, Config*) except +

		World * world
		Config * config

		void run()

 

cdef class pySimulation:
	"""
	The pySimulation class is the main object of revreaddy.
	"""
	cdef Simulation *thisptr
	cdef World *world
	cdef Config *config
	def __cinit__(self):
		self.world = new World();
		self.config = new Config(self.world);
		self.thisptr = new Simulation(self.world, self.config)
	def __dealloc__(self):
		#del self.config
		#del self.world
		del self.thisptr

	property maxTime:
		def __get__(self): return self.config.maxTime
		def __set__(self, maxTime): self.config.maxTime = maxTime
	property temperature:
		def __get__(self): return self.config.temperature
		def __set__(self, temperature): self.config.temperature = temperature
	property timestep:
		def __get__(self): return self.config.timestep
		def __set__(self, timestep): self.config.timestep = timestep
	property cumulativeRuntime:
		def __get__(self): return self.world.cumulativeRuntime
		# there should not be a setter for cumulativeRuntime
		# but here it is
		def __set__(self, cumulativeRuntime):
			self.world.cumulativeRuntime = cumulativeRuntime
	property isPeriodic:
		def __get__(self): return self.config.isPeriodic
		def __set__(self, isPeriodic): self.config.isPeriodic = isPeriodic
	property boxsize:
		def __get__(self): return self.config.boxsize
		def __set__(self, boxsize): self.config.boxsize = boxsize
	property acceptionsDynamics:
		def __get__(self): return self.world.acceptionsDynamics
		def __set__(self, acceptionsDynamics):
			self.world.acceptionsDynamics = acceptionsDynamics
	property rejectionsDynamics:
		def __get__(self): return self.world.rejectionsDynamics
		def __set__(self, rejectionsDynamics):
			self.world.rejectionsDynamics = rejectionsDynamics
	property isReversibleDynamics:
		def __get__(self): return self.config.isReversibleDynamics
		def __set__(self,isReversibleDynamics):
			self.config.isReversibleDynamics = isReversibleDynamics
	property useNeighborList:
		def __get__(self): return self.config.useNeighborList
		def __set__(self,useNeighborList): 
			self.config.useNeighborList=useNeighborList
	property reactionPropagation:
		def __get__(self): return self.config.reactionPropagation
		def __set__(self, reactionPropagation):
			self.config.reactionPropagation = reactionPropagation

	def run(self):
		"""Run the simulation."""
		self.thisptr.run()
	def addParticle(
		self,
		initialPosition=[0.,0.,0.],
		particleTypeId=0):
		self.world.addParticle(
			initialPosition,
			particleTypeId)
	def removeParticle(self, index):
		self.world.removeParticle(index)
	def getPosition(self, index): 
		return self.world.getPosition(index)
	def setPosition(self, index, newPos):
		self.world.setPosition(index, newPos)
	def getTypeId(self, index):
		return self.world.getTypeId(index)
	def setTypeId(self, index, typeId):
		self.world.setTypeId(index, typeId)
	def deleteAllParticles(self):
		self.world.deleteAllParticles()
	def new_Type(
		self,
		name,
		radius,
		diffusionConstant,
		reactionRadius):
		self.config.new_Type(
			name,
			radius,
			diffusionConstant,
			reactionRadius)
	def getNumberOfTypes(self):
		return self.config.getNumberOfTypes()
	def getDictName(self, i):
		return self.config.getDictName(i)
	def getDictRadius(self, i):
		return self.config.getDictRadius(i)
	def getDictDiffusionConstant(self, i):
		return self.config.getDictDiffusionConstant(i)
	def getDictReactionRadius(self, i):
		return self.config.getDictReactionRadius(i)
	def getParticleNumber(self): 
		return self.config.getParticleNumber()
	def writeAllObservablesToFile(self): 
		self.config.writeAllObservablesToFile()
	def writeLastObservableToFile(self):
		self.config.writeLastObservableToFile()
	def showObservables(self):
		return self.config.showObservables()
	def deleteAllObservables(self): 
		self.config.deleteAllObservables()
	def deleteLastObservable(self):
		self.config.deleteLastObservable()
	def new_Trajectory(self, recPeriod, filename): 
		self.config.new_Trajectory(recPeriod, filename)
	# TODO: check sorting of considered along second axis,
	# right: [0,1], [2,4], [1,1]
	# wrong: [1,0], [4,2]
	def new_RadialDistribution(
			self, recPeriod, filename,
			ranges, considered=[[0,0]]):
		self.config.new_RadialDistribution(
			recPeriod, filename,
			ranges, considered)
	def new_MeanSquaredDisplacement(self, recPeriod, filename, particleTypeId):
		self.config.new_MeanSquaredDisplacement(
			recPeriod, filename, particleTypeId)
	def new_ProbabilityDensity(
			self, recPeriod, filename, 
			pTypeId, ranges, coord):
		self.config.new_ProbabilityDensity(
			recPeriod, filename,
			pTypeId, ranges, coord)
	def new_Energy(self, recPeriod, filename):
		self.config.new_Energy(recPeriod, filename)
	def new_Acceptance(self, recPeriod, filename):
		self.config.new_Acceptance(recPeriod, filename)
	def deleteAllGeometries(self):
		self.config.deleteAllGeometries()
	# TODO: check sorting of particleTypeIds, before calling wall constructor
	def new_Wall(self, normal, point, strength, particleTypeIds):
		self.config.new_Wall(normal, point, strength, particleTypeIds)
	def new_DoubleWellZ(self, distanceMinima, strength, particleTypeIds):
		self.config.new_DoubleWellZ(distanceMinima, strength, particleTypeIds)
	def deleteAllForces(self):
		self.config.deleteAllForces()
	def new_SoftRepulsion(self, name, affectedTuple, repulsionStrength):
		self.config.new_SoftRepulsion(name, affectedTuple, repulsionStrength)
	def new_LennardJones(self, name, affectedTuple, epsilon):
		self.config.new_LennardJones(name, affectedTuple, epsilon)
	def getNumberForces(self):
		return self.config.getNumberForces()
	def getForceName(self, index):
		return self.config.getForceName(index)
	def getForceType(self, index):
		return self.config.getForceType(index)
	def getForceAffectedTuple(self, index):
		return self.config.getForceAffectedTuple(index)
	def getForceParameters(self, index):
		return self.config.getForceParameters(index)
		
	# DERIVED FUNCTIONS
	def acceptanceRateDynamics(self):
		if (self.acceptionsDynamics == 0):
			return 0.
		else:
			acc = 1./(1.+ float(self.rejectionsDynamics) \
				/float(self.acceptionsDynamics) )
			return round(acc, 5)

	def efficiency(self):
		return self.acceptanceRate() * self.timestep

	def showTypes(self):
		numberOfTypes = self.getNumberOfTypes()
		print "Number of types:", numberOfTypes
		form = "{:<5}{:<15}{:<10}{:<14}{:<13}"
		print form.format(
			*["Id","Name","Radius","Diffusion-","Reaction-"])
		print form.format(
			*["","","","Constant","Radius"])
		print (5+15+10+14+13)*"-"
		for i in range(numberOfTypes):
			linestr = map(
				str,
				[i, self.getDictName(i), self.getDictRadius(i), 
				 self.getDictDiffusionConstant(i), self.getDictReactionRadius(i) ] )
			print form.format(*linestr)
		return
	
	def showForces(self):
		numberOfForces = self.getNumberForces()
		names = [self.getForceName(i) for i in range(numberOfForces)]
		types = [self.getForceType(i) for i in range(numberOfForces)]
		affectedTuples = [
			self.getForceAffectedTuple(i) for i in range(numberOfForces)
		]
		parameters = [
			self.getForceParameters(i) for i in range(numberOfForces)
		]
		print "Number of forces:", numberOfForces
		form = "{:<15}{:<20}{:<15}{:<15}"
		print form.format(
			*["Name","Type","Affected","Parameters"])
		print (15+20+15+15)*"-"
		for i in range(numberOfForces):
			linestr = map(
				str,
				[names[i], types[i],
				 affectedTuples[i], parameters[i] ] )
			print form.format(*linestr)
		return
