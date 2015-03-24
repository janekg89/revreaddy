import numpy as np
import sim
import h5py
# The state of a simulation consists of four parts:
#   * type dictionary, telling which particle type has which properties
#   * positions, xyz format
#   * reactions, a list of possible reactions
#	* general properties; temperature, boxsize, reversible,
#	  periodic, timestep
# These are put into hdf5 files!

def saveSimulation(filename, simulation):
	fileHandle = h5py.File(filename, "w")

	# retrieve typedict info from simulation
	names = simulation.getDictNames()
	radii = simulation.getDictRadii()
	diffs = simulation.getDictDiffusionConstants()
	reactionRadii = simulation.getDictReactionRadii()
	numberOfTypes = len(names)
	# construct typedict group
	typedict = fileHandle.create_group("typedict")
	typedict.create_dataset("numberOfTypes", (1,), dtype=np.uint64)
	typedict["numberOfTypes"][0] = numberOfTypes
	typedict.create_dataset("names", (numberOfTypes,), dtype="S30")
	typedict.create_dataset("radii", (numberOfTypes,), dtype=np.float64)
	typedict.create_dataset("diffusionConstants", (numberOfTypes,), dtype=np.float64)
	typedict.create_dataset("reactionRadii", (numberOfTypes,), dtype=np.float64)
	for i in range(numberOfTypes):
		typedict["names"][i] = names[i]
		typedict["radii"][i] = radii[i]
		typedict["diffusionConstants"][i] = diffs[i]
		typedict["reactionRadii"][i] = reactionRadii[i]

	# construct xyz group
	xyz = fileHandle.create_group("xyz")
	numberOfParticles = simulation.getParticleNumber()
	xyz.create_dataset("numberOfParticles", (1,), dtype=np.uint64)
	xyz["numberOfParticles"][0] = numberOfParticles
	xyz.create_dataset("typeIds", (numberOfParticles,), dtype=np.uint64)
	xyz.create_dataset("positions", (numberOfParticles,3), dtype=np.float64)
	for i in range(numberOfParticles):
		xyz["typeIds"][i] = simulation.getTypeId(i)
		particlePosition = simulation.getPosition(i)
		xyz["positions"][i,0] = particlePosition[0]
		xyz["positions"][i,1] = particlePosition[1]
		xyz["positions"][i,2] = particlePosition[2]

	# construct properties group
	properties = fileHandle.create_group("properties")
	properties.create_dataset("temperature", (1,), dtype=np.float64)
	properties.create_dataset("boxsize", (1,), dtype=np.float64)
	properties.create_dataset("isReversible", (1,), dtype=np.bool)
	properties.create_dataset("isPeriodic", (1,), dtype=np.bool)
	properties.create_dataset("timestep", (1,), dtype=np.float64)
	properties.create_dataset("repulsionStrength", (1,), dtype=np.float64)
	properties["temperature"][0] = simulation.temperature
	properties["boxsize"][0] = simulation.boxsize
	properties["isReversible"][0] = simulation.isReversible
	properties["isPeriodic"][0] = simulation.isPeriodic
	properties["timestep"][0] = simulation.timestep
	properties["repulsionStrength"][0] = simulation.repulsionStrength

	# construct reactions group
	# ...
	fileHandle.close()

# TODO: Apply hdf5 !
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
