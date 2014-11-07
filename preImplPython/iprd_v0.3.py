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
		consideredParticles = self
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
	
class SimulationType:
	def run(self):
		pass
	def timeLoop(self, activeParticles, maxTime):
		for t in range(0, maxTime):
			consideredParticles = activeParticles.constructConsideredParticles()
			consideredParticles.calculateForces()
			for particle in consideredParticles:
				forceTerm = TIMESTEP * particle.diffConst * particle.force / ( K_BOLTZMANN * TEMPERATURE)
				randomTerm = np.sqrt( 2. * TIMESTEP * particle.diffConst ) * randomNormal()
				particle.move(forceTerm + randomTerm)
				particle.force = np.array([0.,0.,0.])

# The aim of single particle diffusion is to verify the diffusion law:
# 	< ( r(t) - r(0) )**2 > = 6Dt
class SingleParticleDiffusion(SimulationType):
	def run(self, maxRepetitions, maxTime, maxCounts):
		endDistances = np.zeros((maxRepetitions, maxCounts))
		measurementTimes = np.arange(0.,float(maxTime)*TIMESTEP, float(maxTime) * TIMESTEP / float(maxCounts))
		theParticle = Particle("the diffusing particle")
		activeParticles = ActiveParticles()
		activeParticles.extend([theParticle])		
		for repetition in range(0, maxRepetitions):
			theParticle.position = np.array([0.,0.,0.])
			startPosition = cp.deepcopy(theParticle.position)
			for counter in range(0, maxCounts):
				endDistances[repetition, counter] = squaredDistance(startPosition, theParticle.position)
				self.timeLoop(activeParticles, maxTime/maxCounts)
		return (endDistances, measurementTimes)
				
def randomNormal(): return np.random.normal(0., 1., 3)

def squaredDistance(r1, r2):
	return ( r2[0] - r1[0] )**2 + ( r2[1] - r1[1] )**2 + ( r2[2] - r1[2] )**2


"""---------------------------------------------------------------------------------"""
"""---------------------------------------------------------------------------------"""
"""---------------------------------------------------------------------------------"""

TIMESTEP 		= 1.
K_BOLTZMANN 	= 1.
TEMPERATURE 	= 1.

simulation = SingleParticleDiffusion()
(endDistances, measurementTimes) = simulation.run(1000, 200, 20)
print endDistances

distances = np.zeros(len(endDistances[0]))

for time in range(0,len(endDistances[0])):
	meanDistance = 0.
	for rep in range(0, len(endDistances)):
		meanDistance += endDistances[rep,time]
	distances[time] = meanDistance / float(len(endDistances))

plt.plot(measurementTimes, distances)
plt.show()
