import numpy as np
import matplotlib.pyplot as plt

# activeParticles keeps track of which particles are 
# involved in the simlutaion. For these particles
# forces and reactions will be considered and time is
# integrated. 
global activeParticles
activeParticles = []

def shiftActiveIndicesFrom(index):
	for i in range(index, len(activeParticles)):
		activeParticles[i].activeIndex -= 1

class Particle:
	def __init__(self, name):
		self.name 		= name
		self.position 	= np.array([0.,0.,0.])
		self.timestep 	= 0
		self.activeIndex= 0
		self.radius		= 1.
		self.type		= 'misc'
	
	def move(self, delta):
		self.position += delta
	
	def appear(self):
		activeParticles.append(self)
		self.activeIndex = len(activeParticles) - 1
	
	def disappear(self):
		del activeParticles[self.activeIndex]
		shiftActiveIndicesFrom(self.activeIndex)
	
	def getName(self):
		return self.name
	
	def getActiveIndex(self):
		return self.activeIndex

a = Particle('h20')
b = Particle('co2')
c = Particle('h')
d = Particle('c')

a.appear()
b.appear()
c.appear()
d.appear()

print map(Particle.getActiveIndex, activeParticles)
print map(Particle.getName, activeParticles)

b.disappear()

print map(Particle.getActiveIndex, activeParticles)
print map(Particle.getName, activeParticles)

print b.name
