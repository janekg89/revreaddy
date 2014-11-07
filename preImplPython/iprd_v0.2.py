import numpy as np
import matplotlib.pyplot as plt
import copy as cp

class Particle:
	def __init__(self, name):
		self.name 		= name
		self.position 	= np.array([0.,0.,0.])
		self.timestep 	= 0
		self.count		= 0
		self.skip		= 0
		self.radius		= 1.
		self.diffConst	= 1.
		self.force		= np.array([0.,0.,0.])
		self.type		= 'misc'
	
	def move(self, delta):
		self.position += delta
	
	def getName(self): return self.name

class ActiveParticles(list):
	def constructConsideredParticles(self):
		self.findTimeskips()
		self.orderByCount()
		# this has to become a more efficient function in C++.
		consideredParticles = cp.deepcopy(self)
		# consideredParticles = [particle for particle in self if particle.count == 0]
		return consideredParticles
	
	# findTimeskips will iterate over ActiveParticles to 
	# set the skippable timesteps for every particle.
	# Also will set the count initially equal to the timeskips.
	def findTimeskips(self):
		pass
	
	# orderByCount will order the activeParticles list according
	# to the particles' count value.	
	def orderByCount(self):
		pass

	def calculateForces(self):
		for i in range(0,len(self)):
			for j in range(i,len(self)):
				#currentForce = softcoreForce(self[i], self[j])
				currentForce = np.array([0.,0.,0.])
				self[i].force += currentForce
				self[j].force -= currentForce
	

def randomNormal(): return np.random.normal(0., 1., 3)

"""---------------------------------------------------------------------------------"""
"""---------------------------------------------------------------------------------"""
"""---------------------------------------------------------------------------------"""

TIMESTEP 	= 1.
K_BOLTZMANN = 1.
TEMPERATURE = 1.

activeParticles = ActiveParticles()
consideredParticles = ActiveParticles()

particleA = Particle('particleA')
particleB = Particle('particleB')

activeParticles.extend([particleA])
activeParticles.extend([particleB])
#print map(Particle.getName, activeParticles)

print activeParticles[0].position
activeParticles[0].move(np.array([1.,2.,3.]))
print activeParticles[0].position
print activeParticles[1].position
print particleA.position

"""
for t in range(0,10):
	consideredParticles = activeParticles.constructConsideredParticles()
	consideredParticles.calculateForces()
	for particle in consideredParticles:
		forceTerm = TIMESTEP * particle.diffConst * particle.force / ( K_BOLTZMANN * TEMPERATURE)
		randomTerm = np.sqrt( 2. * TIMESTEP * particle.diffConst ) * randomNormal()
		particle.move(forceTerm + randomTerm)
		particle.force = np.array([0.,0.,0.])
			

"""
