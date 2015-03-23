import numpy as np
import sim
# TODO: apply the following idea!
# The state of a simulation consists of four parts:
#   * type dictionary .dic, telling which particle type has which properties
#   * positions, a .xyz file containing a single frame
#   * reactions .rea, a list of possible reactions
#	* general properties .pro, temperature, boxsize, reversible,
#	  periodic, timestep

def saveSimulation(filename, simulation):
	# construct .dic file
	names = simulation.getDictNames()
	radii = simulation.getDictRadii()
	diffs = simulation.getDictDiffusionConstants()
	reactionRadii = simulation.getDictReactionRadii()
	numberOfTypes = len(names)
	dicFile = open(filename+".dic", "w")
	dicFile.write("Id\tName\tRadius\tDiffusionConstant\tReactionRadius\n")
	for i in range(numberOfTypes):
		line  = str(i)
		line += str(names[i]) + "\t" + str(radii[i]) + "\t"
		line += str(diffs[i]) + "\t" + str(reactionRadii[i]) + "\n"
		dicFile.write(line)
	dicFile.close()

	# construct .xyz file
	xyzFile = open(filename+".xyz", "w")
	numberOfParticles = simulation.getParticleNumber()
	xyzFile.write(str(numberOfParticles) + "\n")
	xyzFile.write("#particleTypeId\tx\ty\tz\n") # comment line
	for i in range(numberOfParticles):
		particleType = simulation.getTypeId(i)
		particlePosition = simulation.getPosition(i)
		line  = str(particleType) + "\t"
		line += str(particlePosition[0]) + "\t"
		line += str(particlePosition[1]) + "\t"
		line += str(particlePosition[2]) + "\n"
	xyzFile.close()

	# construct .pro file
	proFile = open(filename+".pro", "w")
	line  = "Temperature\tBoxsize\tisReversible\t"
	line += "isPeriodic\tTimestep\tRepulsionStrength\n"
	proFile.write(line)
	line  = str(simulation.temperature) + "\t"
	line += str(simulation.boxsize) + "\t"
	# bools are converted to ints (0 or 1) because the
	# function bool() returns True when given the string "False"
	# since "False" is a non-zero expression (thinking of bits)
	line += str(int(simulation.isReversible)) + "\t"
	line += str(int(simulation.isPeriodic)) + "\t"
	line += str(simulation.timestep) + "\t"
	line += str(simulation.repulsionStrength) + "\n"
	proFile.write(line)
	proFile.close()

	# construct .rea file 
	# ...

# TODO: apply idea with .dic, .xyz, .pro and .rea files
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
