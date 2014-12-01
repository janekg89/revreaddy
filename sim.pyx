# distutils: language = c++
# distutils: extra_compile_args = -std=c++11
# distutils: sources = src/Simulation.cpp src/Particle.cpp src/Random.cpp src/Potential.cpp src/Filehandler.cpp
# distutils: include_dirs = /usr/include /usr/local/include include

from libcpp.vector cimport vector
cimport numpy as np

cdef extern from "Simulation.h":
	cdef cppclass Simulation:
		unsigned long int maxTime
		unsigned long int getMaxTime()
		void addParticle()
		vector[double] getParticlePosition()
		void run()

cdef class pySimulation:
	cdef Simulation *thisptr
	def __cinit__(self):
		self.thisptr = new Simulation()
	def __dealloc__(self):
		del self.thisptr
	def getMaxTime(self):
		return self.thisptr.getMaxTime()
	def addParticle(self):
		self.thisptr.addParticle()
	def getParticlePosition(self):
		return self.thisptr.getParticlePosition()
	def run(self):
		self.thisptr.run()
	property maxTime:
		def __get__(self): return self.thisptr.maxTime
		def __set__(self, maxTime): self.thisptr.maxTime = maxTime
