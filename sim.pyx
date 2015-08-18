#!python
#cython: embedsignature=True
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp cimport bool
cimport numpy as np
import time

cdef extern from "Simulation.h":
	cdef cppclass Simulation:
		Simulation() except +
		unsigned long int maxTime
		double temperature
		double timestep
		double cumulativeRuntime
		bool isPeriodic
		double boxsize
		unsigned long int acceptionsDynamics
		unsigned long int rejectionsDynamics
		bool isReversibleDynamics
		bool isReversibleReactions
		bool useNeighborList

		void run()
		void addParticle(vector[double], unsigned int)
		vector[double] getPosition(int)
		void           setPosition(int, vector[double])
		unsigned int getTypeId(int)
		void         setTypeId(int, unsigned int)
		void new_Type(string, double, double, double)
		vector[string] getDictNames()
		vector[double] getDictRadii()
		vector[double] getDictDiffusionConstants()
		vector[double] getDictReactionRadii()
		int getParticleNumber()
		void deleteAllParticles()
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
		void new_MeanSquaredDisplacement(unsigned long int, string, unsigned int)
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

cdef class pySimulation:
	"""
	The pySimulation class is the main object of revreaddy.
	"""
	cdef Simulation *thisptr
	def __cinit__(self):
		self.thisptr = new Simulation()
	def __dealloc__(self):
		del self.thisptr
	def addParticle(
		self,
		initialPosition=[0.,0.,0.],
		particleTypeId=0):
		self.thisptr.addParticle(
			initialPosition,
			particleTypeId)
	def run(self):
		"""Run the simulation."""
		t1 = time.clock()
		self.thisptr.run()
		t2 = time.clock()
		print "Time needed:", t2 - t1, "seconds"
	def getPosition(self, index): 
		return self.thisptr.getPosition(index)
	def setPosition(self, index, newPos):
		self.thisptr.setPosition(index, newPos)
	def getTypeId(self, index):
		return self.thisptr.getTypeId(index)
	def setTypeId(self, index, typeId):
		self.thisptr.setTypeId(index, typeId)
	def new_Type(
		self,
		name,
		radius,
		diffusionConstant,
		reactionRadius):
		self.thisptr.new_Type(
			name,
			radius,
			diffusionConstant,
			reactionRadius)
	def getDictNames(self):
		return self.thisptr.getDictNames()
	def getDictRadii(self):
		return self.thisptr.getDictRadii()
	def getDictDiffusionConstants(self):
		return self.thisptr.getDictDiffusionConstants()
	def getDictReactionRadii(self):
		return self.thisptr.getDictReactionRadii()
	def getParticleNumber(self): 
		return self.thisptr.getParticleNumber()
	def deleteAllParticles(self):
		self.thisptr.deleteAllParticles()
	def writeAllObservablesToFile(self): 
		self.thisptr.writeAllObservablesToFile()
	def writeLastObservableToFile(self):
		self.thisptr.writeLastObservableToFile()
	def showObservables(self):
		return self.thisptr.showObservables()
	def deleteAllObservables(self): 
		self.thisptr.deleteAllObservables()
	def deleteLastObservable(self):
		self.thisptr.deleteLastObservable()
	def new_Trajectory(self, recPeriod, filename): 
		self.thisptr.new_Trajectory(recPeriod, filename)
	# TODO: check sorting of considered along second axis,
	# right: [0,1], [2,4], [1,1]
	# wrong: [1,0], [4,2]
	def new_RadialDistribution(
			self, recPeriod, filename,
			ranges, considered=[[0,0]]):
		self.thisptr.new_RadialDistribution(
			recPeriod, filename,
			ranges, considered)
	def new_MeanSquaredDisplacement(self, recPeriod, filename, particleTypeId):
		self.thisptr.new_MeanSquaredDisplacement(
			recPeriod, filename, particleTypeId)
	def new_ProbabilityDensity(
			self, recPeriod, filename, 
			pTypeId, ranges, coord):
		self.thisptr.new_ProbabilityDensity(
			recPeriod, filename,
			pTypeId, ranges, coord)
	def new_Energy(self, recPeriod, filename):
		self.thisptr.new_Energy(recPeriod, filename)
	def new_Acceptance(self, recPeriod, filename):
		self.thisptr.new_Acceptance(recPeriod, filename)
	def deleteAllGeometries(self):
		self.thisptr.deleteAllGeometries()
	# TODO: check sorting of particleTypeIds, before calling wall constructor
	def new_Wall(self, normal, point, strength, particleTypeIds):
		self.thisptr.new_Wall(normal, point, strength, particleTypeIds)
	def new_DoubleWellZ(self, distanceMinima, strength, particleTypeIds):
		self.thisptr.new_DoubleWellZ(distanceMinima, strength, particleTypeIds)
	def deleteAllForces(self):
		self.thisptr.deleteAllForces()
	def new_SoftRepulsion(self, name, affectedTuple, repulsionStrength):
		self.thisptr.new_SoftRepulsion(name, affectedTuple, repulsionStrength)
	def new_LennardJones(self, name, affectedTuple, epsilon):
		self.thisptr.new_LennardJones(name, affectedTuple, epsilon)
	def getNumberForces(self):
		return self.thisptr.getNumberForces()
	def getForceName(self, index):
		return self.thisptr.getForceName(index)
	def getForceType(self, index):
		return self.thisptr.getForceType(index)
	def getForceAffectedTuple(self, index):
		return self.thisptr.getForceAffectedTuple(index)
	def getForceParameters(self, index):
		return self.thisptr.getForceParameters(index)
		
	property maxTime:
		def __get__(self): return self.thisptr.maxTime
		def __set__(self, maxTime): self.thisptr.maxTime = maxTime
	property temperature:
		def __get__(self): return self.thisptr.temperature
		def __set__(self, temperature): self.thisptr.temperature = temperature
	property timestep:
		def __get__(self): return self.thisptr.timestep
		def __set__(self, timestep): self.thisptr.timestep = timestep
	property cumulativeRuntime:
		def __get__(self): return self.thisptr.cumulativeRuntime
		# there should not be a setter for cumulativeRuntime
		# but here it is
		def __set__(self, cumulativeRuntime):
			self.thisptr.cumulativeRuntime = cumulativeRuntime
	property isPeriodic:
		def __get__(self): return self.thisptr.isPeriodic
		def __set__(self, isPeriodic): self.thisptr.isPeriodic = isPeriodic
	property boxsize:
		def __get__(self): return self.thisptr.boxsize
		def __set__(self, boxsize): self.thisptr.boxsize = boxsize
	property acceptionsDynamics:
		def __get__(self): return self.thisptr.acceptionsDynamics
		def __set__(self, acceptionsDynamics):
			self.thisptr.acceptionsDynamics = acceptionsDynamics
	property rejectionsDynamics:
		def __get__(self): return self.thisptr.rejectionsDynamics
		def __set__(self, rejectionsDynamics):
			self.thisptr.rejectionsDynamics = rejectionsDynamics
	property isReversibleDynamics:
		def __get__(self): return self.thisptr.isReversibleDynamics
		def __set__(self,isReversibleDynamics):
			self.thisptr.isReversibleDynamics = isReversibleDynamics
	property useNeighborList:
		def __get__(self): return self.thisptr.useNeighborList
		def __set__(self,useNeighborList): 
			self.thisptr.useNeighborList=useNeighborList

	# derived functions
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
		names = self.getDictNames()
		radii = self.getDictRadii()
		diffs = self.getDictDiffusionConstants()
		reactionRadii = self.getDictReactionRadii()
		numberOfTypes = len(names)
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
				[i, names[i], radii[i], 
				 diffs[i], reactionRadii[i] ] )
			print form.format(*linestr)
		return
	
	def showForces(self):
		numberOfForces = self.getNumberForces()
		names = [self.getForceName(i) for i in range(numberOfForces)]
		types = [self.getForceType(i) for i in range(numberOfForces)]
		affectedTuples = [self.getForceAffectedTuple(i) for i in range(numberOfForces)]
		parameters = [self.getForceParameters(i) for i in range(numberOfForces)]
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
