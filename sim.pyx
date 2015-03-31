# distutils: language = c++
# distutils: extra_compile_args = -std=c++11 
# distutils: sources = src/Observable.cpp src/Simulation.cpp src/Particle.cpp src/Random.cpp src/Force.cpp src/Trajectory.cpp src/TrajectorySingle.cpp src/RadialDistribution.cpp src/utils.cpp src/MeanSquaredDisplacement.cpp src/ProbabilityDensity.cpp src/TypeDict.cpp
# distutils: include_dirs = /usr/include /usr/local/include include
# distutils: libraries = m gsl gslcblas

from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp cimport bool
cimport numpy as np

cdef extern from "Simulation.h":
	cdef cppclass Simulation:
		Simulation(bool) except +
		unsigned long int maxTime
		double temperature
		double timestep
		double cumulativeRuntime
		bool isPeriodic
		double boxsize
		double repulsionStrength
		bool verbose
		unsigned long int acceptions
		unsigned long int rejections
		bool isReversible

		void addParticle(vector[double], unsigned int)
		void run()
		vector[double] getPosition(int)
		void           setPosition(int, vector[double])
		unsigned int getTypeId(int)
		void         setTypeId(int, unsigned int)
		void new_Type(string, double, double, double, unsigned int)
		vector[string] getDictNames()
		vector[double] getDictRadii()
		vector[double] getDictDiffusionConstants()
		vector[double] getDictReactionRadii()
		vector[uint]   getDictForceTypes()
		int getParticleNumber()
		void deleteAllParticles()
		void writeAllObservablesToFile()
		string showObservables()
		void deleteAllObservables()
		void new_Trajectory(string)
		void new_TrajectorySingle()
		vector[vector[double]] getTrajectorySingle()
		void new_RadialDistribution(string,vector[double],vector[vector[uint]])
		void new_MeanSquaredDisplacement(string, unsigned int)
		void new_ProbabilityDensity(
			string, unsigned int, 
			vector[double], unsigned int)

cdef class pySimulation:
	cdef Simulation *thisptr
	def __cinit__(self, hasDefaultTypes=True):
		self.thisptr = new Simulation(hasDefaultTypes)
	def __dealloc__(self):
		del self.thisptr
	def addParticle(
		self,
		initialPosition=[0.,0.,0.],
		particleTypeId=2):
		self.thisptr.addParticle(
			initialPosition,
			particleTypeId)
	def run(self):
		self.thisptr.run()
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
		reactionRadius,
		forceType):
		self.thisptr.new_Type(
			name,
			radius,
			diffusionConstant,
			reactionRadius,
			forceType)
	def getDictNames(self):
		return self.thisptr.getDictNames()
	def getDictRadii(self):
		return self.thisptr.getDictRadii()
	def getDictDiffusionConstants(self):
		return self.thisptr.getDictDiffusionConstants()
	def getDictReactionRadii(self):
		return self.thisptr.getDictReactionRadii()
	def getDictForceTypes(self):
		return self.thisptr.getDictForceTypes()
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
	def new_RadialDistribution(self, filename, ranges, considered=[[2,2]]):
		self.thisptr.new_RadialDistribution(filename, ranges, considered)
	def new_MeanSquaredDisplacement(self, filename, particleTypeId):
		self.thisptr.new_MeanSquaredDisplacement(filename, particleTypeId)
	def new_ProbabilityDensity(self, filename, pTypeId, ranges, coord):
		self.thisptr.new_ProbabilityDensity(filename, pTypeId, ranges, coord)
		
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
	property isReversible:
		def __get__(self): return self.thisptr.isReversible
		def __set__(self,isReversible): self.thisptr.isReversible=isReversible

	# derived functions
	def acceptanceRate(self):
		if (self.acceptions == 0):
			return 0.
		else:
			acc = 1./(1.+ float(self.rejections)/float(self.acceptions) )
			return round(acc, 5)
	def efficiency(self):
		return self.acceptanceRate() * self.timestep

	def showTypes(self):
		names = self.getDictNames()
		radii = self.getDictRadii()
		diffs = self.getDictDiffusionConstants()
		reactionRadii = self.getDictReactionRadii()
		forceTypes = self.getDictForceTypes()
		numberOfTypes = len(names)
		print "Number of types:", numberOfTypes
		form = "{:<5}{:<15}{:<10}{:<14}{:<13}{:<8}"
		print form.format(
			*["Id","Name","Radius","Diffusion-","Reaction-","Force-"])
		print form.format(
			*["","","","Constant","Radius","Type"])
		print (5+15+10+14+13+8)*"-"
		for i in range(numberOfTypes):
			linestr = map(
				str,
				[i, names[i], radii[i], 
				 diffs[i], reactionRadii[i], forceTypes[i] ] )
			print form.format(*linestr)
		
# this dict() sets the type ids used in the program
# the first three keys and values (none,0),(lj,1),(soft,2)
# should NOT be changed.
typeStringToId = dict()
typeStringToId["none"] = 0
typeStringToId["lj"]   = 1
typeStringToId["soft"] = 2
