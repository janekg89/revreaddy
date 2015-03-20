import numpy as np
import sim
# TODO: apply the following idea!
# The state of a simulation consists of two parts:
#   * type dictionary .dic, telling which particle type has which properties
#   * positions, a .xyz file containing a single frame
#   * reactions .rea, a list of possible reactions

def saveSimulation(filename, simulation):
	# system contains [0] particlenumber, [1] boxsize, [2] temperature
	system = np.zeros(3)
	system[0] = float(simulation.getParticleNumber())
	system[1] = simulation.boxsize
	system[2] = simulation.temperature
	system    = map(str, system)
	system    = map(lambda x: x+"\t", system)
	system    = "".join(system)
	system    = system[:-1]
	# particles contain per row:
	# [0] typeId, [1] radius, [2] diffusionConstant,,
	# [3] x, [4] y, [5] z - coordinates
	N = simulation.getParticleNumber()
	particles = np.zeros((N, 6))
	for i in range(N):
		particles[i,0] = float(simulation.getTypeId(i))
		particles[i,1] = simulation.getRadius(i)
		particles[i,2] = simulation.getDiffusionConstant(i)
		particles[i,3] = simulation.getPosition(i)[0]
		particles[i,4] = simulation.getPosition(i)[1]
		particles[i,5] = simulation.getPosition(i)[2]
	np.savetxt(filename, particles,  header=system, comments="")

	
def loadSimulation(filename):
	with open(filename, 'r') as fHandle:
		firstLine  =  fHandle.readline()
	system = np.fromstring(firstLine, sep="\t")
	simulation             = sim.pySimulation()
	simulation.boxsize     = system[1]
	simulation.temperature = system[2]

	N = int(system[0])
	particles = np.loadtxt(filename, skiprows=1)
	for i in range(N):
		simulation.addParticle(
			initialPosition   = [
				particles[i,3],
				particles[i,4],
				particles[i,5]
			],
			particleTypeId    = int(particles[i,0]),
			radius            = particles[i,1],
			diffusionConstant = particles[i,2]
		)
	return simulation
