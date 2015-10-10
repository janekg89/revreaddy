#!python
#cython: embedsignature=True
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp cimport bool
cimport numpy as np
import numpy as np
import scipy.integrate
import time
import logging

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
		unsigned numberBoxes
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
		void new_Acceptance(unsigned long int, string, bool)
		void new_ParticleNumbers(unsigned long, string, unsigned)
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
		double getForceCutoff(unsigned)
		void deleteAllReactions()
		void new_Conversion(string, unsigned int, unsigned int, double, double)
		void new_Fusion(string, unsigned, unsigned, unsigned, double, double)
		void configure_Fusion(
			unsigned, vector[uint], double, 
			double, double, double, double, double)
		unsigned getNumberReactions()
		string getReactionName(unsigned)
		string getReactionType(unsigned)
		vector[uint] getReactionForwardTypes(unsigned)
		vector[uint] getReactionBackwardTypes(unsigned)
		double getReactionForwardRate(unsigned)
		double getReactionBackwardRate(unsigned)

	cdef cppclass Simulation:
		Simulation() except +

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
		# Simulation allocates World and Config for itself
		self.thisptr = new Simulation()
		self.world = self.thisptr.world
		self.config = self.thisptr.config
		logging.info("pySimulation constructed.")
	def __dealloc__(self):
		del self.thisptr
		logging.info("pySimulation released.")

	property maxTime:
		def __get__(self): return self.config.maxTime
		def __set__(self, maxTime): self.config.maxTime = maxTime
	property temperature:
		def __get__(self): return self.config.temperature
		def __set__(self, temperature): self.config.temperature = temperature
	property kBoltzmann:
		def __get__(self): return self.config.kBoltzmann
		def __set__(self, kBoltzmann): self.config.kBoltzmann = kBoltzmann
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
		# TODO when boxize is varied, check if particles are still inside
		def __set__(self, boxsize): self.config.boxsize = boxsize
	property acceptionsDynamics:
		def __get__(self): return self.world.acceptionsDynamics
		def __set__(self, acceptionsDynamics):
			self.world.acceptionsDynamics = acceptionsDynamics
	property rejectionsDynamics:
		def __get__(self): return self.world.rejectionsDynamics
		def __set__(self, rejectionsDynamics):
			self.world.rejectionsDynamics = rejectionsDynamics
	property acceptionsReactions:
		def __get__(self): return self.world.acceptionsReactions
		def __set__(self, acceptionsReactions):
			self.world.acceptionsReactions = acceptionsReactions
	property rejectionsReactions:
		def __get__(self): return self.world.rejectionsReactions
		def __set__(self, rejectionsReactions):
			self.world.rejectionsReactions = rejectionsReactions
	property isReversibleDynamics:
		def __get__(self): return self.config.isReversibleDynamics
		def __set__(self,isReversibleDynamics):
			self.config.isReversibleDynamics = isReversibleDynamics
	property isReversibleReactions:
		def __get__(self): return self.config.isReversibleReactions
		def __set__(self,isReversibleReactions):
			self.config.isReversibleReactions = isReversibleReactions
	property useNeighborList:
		def __get__(self): return self.config.useNeighborList
		def __set__(self,useNeighborList): 
			self.config.useNeighborList=useNeighborList
	property numberBoxes:
		def __get__(self): return self.config.numberBoxes
		def __set__(self, numberBoxes):
			self.config.numberBoxes = numberBoxes
	property reactionPropagation:
		def __get__(self): return self.config.reactionPropagation
		def __set__(self, reactionPropagation):
			self.config.reactionPropagation = reactionPropagation

	def run(self, maxTime=None, timestep=None):
		"""Run the simulation."""
		self.configureAllReactions()
		self.checkIfNeighborlist()
		if (maxTime != None): self.maxTime = maxTime
		if (timestep != None): self.timestep = timestep
		logging.info("Run() with timestep: " + str(self.timestep) + \
			" and maxTime: " + str(self.maxTime))
		logging.info("Started at " + \
			str(self.world.cumulativeRuntime) + " cumulative runtime.")
		timer = time.clock()
		self.thisptr.run()
		timer = time.clock() - timer
		logging.info("Stopped at " + \
			str(self.world.cumulativeRuntime) + " cumulative runtime.")
		logging.info("Needed " + str(timer) + " s for computation.")
	def addParticle(
		self,
		initialPosition=[0.,0.,0.],
		particleTypeId=0):
		if (particleTypeId >= self.getNumberOfTypes()):
			logging.error(
				"Particle type does not exist. Particle is not created.")
			return
		distances = np.array( map(abs, initialPosition) )
		distances -= self.boxsize/2.
		outsideBox = map(lambda x: x>0, distances)
		outsideBox = reduce(lambda x,y: x or y, outsideBox)
		if ( outsideBox ):
			logging.error(
				"Initial position would be outside of box. " + \
				"Particle is not created.")
			return
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
	def new_Acceptance(self, recPeriod, filename, reactionsOrDynamics):
		self.config.new_Acceptance(recPeriod, filename, reactionsOrDynamics)
	def new_ParticleNumbers(self, recPeriod, filename, particleTypeId):
		self.config.new_ParticleNumbers(recPeriod, filename, particleTypeId)
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
	def getForceCutoff(self, index):
		return self.config.getForceCutoff(index)
	def deleteAllReactions(self):
		self.config.deleteAllReactions()
	def new_Conversion(self, name, forwardType,
		backwardType, forwardRate, backwardRate):
		if (forwardType >= self.getNumberOfTypes()):
			logging.error(
				"ForwardType does not exist. Reaction not created.")
			return
		if (backwardType >= self.getNumberOfTypes()):
			logging.error(
				"BackwardType does not exist. Reaction not created.")
			return
		self.config.new_Conversion(name, forwardType,
			backwardType, forwardRate, backwardRate)
	def new_Fusion(self, name, forwardTypeA, forwardTypeB,
		backwardTypeC, forwardRate, backwardRate):
		self.config.new_Fusion(name, forwardTypeA, forwardTypeB,
			backwardTypeC, forwardRate, backwardRate)
	def configure_Fusion(self, reactionIndex, interactionsIndices,
		inversePartition, maxDistr, radiiSum, reactionRadiiSum,
		meanDistr, inverseTemperature):
		self.config.configure_Fusion(reactionIndex, interactionsIndices,
			inversePartition, maxDistr, radiiSum, reactionRadiiSum,
			meanDistr, inverseTemperature)
	def getNumberReactions(self):
		return self.config.getNumberReactions()
	def getReactionName(self, index):
		return self.config.getReactionName(index)
	def getReactionType(self, index):
		return self.config.getReactionType(index)
	def getReactionForwardTypes(self, index):
		return self.config.getReactionForwardTypes(index)
	def getReactionBackwardTypes(self, index):
		return self.config.getReactionBackwardTypes(index)
	def getReactionForwardRate(self, index):
		return self.config.getReactionForwardRate(index)
	def getReactionBackwardRate(self, index):
		return self.config.getReactionBackwardRate(index)
		
	# DERIVED FUNCTIONS
	def configureAllReactions(self):
		logging.info("Configuring all reactions ...")
		# search all reactions and find affected particleTypes
		# from those find relevant interactions and calculate
		# inverse partition function, maximum of distribution
		# mean of distribution, sum of reaction radii and inverse
		# temperature
		for i in range(self.getNumberReactions()):
			if (self.getReactionType(i) == "Fusion"):
				self.configureFusion(i) 

	def configureFusion(self, reactionIndex):
		# Fusion is A + B <-> C, 
		# hence len(forwardTypes) = 2 and len(backwardTypes) = 1
		# only forwardTypes might have interactions
		logging.info("Configure Fusion reaction #" + str(reactionIndex))
		forwardTypes = self.getReactionForwardTypes(reactionIndex)
		aType = forwardTypes[0]
		bType = forwardTypes[1]
		interactions = [] # indices of the interactions
		for i in range(self.getNumberForces()):
			affectedTuple = self.getForceAffectedTuple(i)
			if ( 
				( (aType == affectedTuple[0])
				and (bType == affectedTuple[1]) )
				or 
				( (aType == affectedTuple[1]) 
				and (bType == affectedTuple[0]) ) ):
				interactions += [i]
		# radii of particles A and B
		radiusA = self.getDictRadius(aType)
		radiusB = self.getDictRadius(bType)
		radiiSum = radiusA + radiusB
		# define the energy function for these 
		# interactions as a sum of the single terms
		def energy(distance):
			result = 0.
			for inter in interactions:
				if ( self.getForceType(inter) == "SoftRepulsion"):
					strength = self.getForceParameters(inter)[0]
					result += self.softRepulsionEnergy(distance, radiiSum, strength)
				elif (self.getForceType(inter) == "LennardJones"):
					epsilon = self.getForceParameters(inter)[0]
					result += self.lennardJonesEnergy(distance, radiiSum, epsilon)
			return result

		# reaction radii of particles A and B
		reactionRadiusA = self.getDictReactionRadius(aType)
		reactionRadiusB = self.getDictReactionRadius(bType)
		reactionRadiiSum = reactionRadiusA + reactionRadiusB
		inverseTemperature = 1. / (self.temperature * self.kBoltzmann)
		integrand = lambda x: \
			x * np.exp( -inverseTemperature * energy( x*reactionRadiiSum ))
		unitRange = np.arange(0., 1.001, 0.001)
		yRange = map(integrand, unitRange)
		partitionFunction = scipy.integrate.simps(yRange, unitRange)
		inversePartition = 1. / partitionFunction
		# define the distribution from which the uniform numbers for
		# particle distances will be drawn
		def distribution(x):
			result = np.exp(-inverseTemperature*energy(x*reactionRadiiSum))
			return x * inversePartition * result

		# now only missing are the distributions' mean value and its maximum
		# in the interval [0,1]
		maxDistr = max( map(distribution, unitRange) )
		integrand = lambda x: x * distribution(x)
		yRange = map(integrand, unitRange)
		meanDistr = scipy.integrate.simps(yRange, unitRange) 

		self.config.configure_Fusion(
			reactionIndex,
			interactions,
			inversePartition,
			maxDistr,
			radiiSum,
			reactionRadiiSum,
			meanDistr,
			inverseTemperature)
		logging.info(
			"Configured Fusion reaction with parameters: "+\
			"reactionIndex " + str(reactionIndex) +\
			", interactions " + " ".join( map(str,interactions) ) +\
			", inversePartition " + str(inversePartition) +\
			", maxDistr " + str(maxDistr) +\
			", radiiSum " + str(radiiSum) +\
			", reactionRadiiSum " + str(reactionRadiiSum) +\
			", meanDistr " + str(meanDistr) +\
			", inverseTemperature " + str(inverseTemperature)
		)

	def softRepulsionEnergy(self, distance, radiiSum, strength):
		if ( distance > radiiSum ): return 0.
		return strength * ( distance - radiiSum )**2


	def lennardJonesEnergy(self, distance, radiiSum, epsilon):
		if ( distance > (2.5*radiiSum) ): return 0.
		# sigma = 2.**(-1/6) * radiiSum
		sigma = 0.8908987181403393 * radiiSum
		return 4.*epsilon*( (sigma/distance)**12 - (sigma/distance)**6 )

	def checkIfNeighborlist(self):
		# We want to find out if neighborlist pays off. 
		# First determine the minimal size of subboxes,
		# therefore check all interaction distances and
		# all reaction radii.	
		minimalLength = self.boxsize
		# check interaction distances
		for i in range(self.getNumberForces()):
			if (self.getForceCutoff(i) < minimalLength):
				minimalLength = self.getForceCutoff(i)
		# check reaction radii combinations (i,j)
		for i in range(self.getNumberOfTypes()):
			for j in range(i, self.getNumberOfTypes()):
				R=self.getDictReactionRadius(i)+self.getDictReactionRadius(j)
				if (R < minimalLength):
					minimalLength = R
		counter, numberBoxes = 1., 0
		while ((self.boxsize / counter) > minimalLength):
			numberBoxes += 1
			counter += 1.
		# if n=3 we will have 9 subboxes of length L/n, which
		# will result in having to check every box. This is
		# as inefficient as double looping. So:
		# ONLY construct neighborlist if we have at least 16
		# subboxes or n>3
		if ((numberBoxes > 3) and self.useNeighborList):
			self.useNeighborList = True # obsolete, but who cares
			self.numberBoxes = numberBoxes
			logging.info("useNeighborList enabled with " + str(numberBoxes) \
				+ " number of boxes along each axis.")
		else:
			self.useNeighborList = False

	# UTILITY
	def acceptanceRateDynamics(self):
		if (self.acceptionsDynamics == 0):
			return 0.
		else:
			acc = 1./(1.+ float(self.rejectionsDynamics) \
				/float(self.acceptionsDynamics) )
			return round(acc, 5)

	def acceptanceRateReactions(self):
		if (self.acceptionsReactions == 0):
			return 0.
		else:
			acc = 1./(1.+ float(self.rejectionsReactions) \
				/float(self.acceptionsReactions) )
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

	def howManyParticles(self, particleType):
		counter = 0
		for i in range(0, self.getParticleNumber()):
			if (self.getTypeId(i) == particleType): counter += 1
		return (str(counter) + "/" + str(self.getParticleNumber()))