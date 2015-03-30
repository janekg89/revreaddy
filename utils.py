"""
The utils module provides necessary tools when working with revreaddy. 
Saving and loading the state of a simulation is one of the tasks of
this module. At a given time during your simulation the state of the 
system can be partitioned into four parts:
	* type dictionary, telling which particle type has which properties
	* positions, a structure holding a typeId and xyz coordinates for
	  every particle
	* reactions, a list of possible reactions
	* general properties of the system

.. note::
	Observables are not part of the underlying system itself 
	and will not be saved.

These four are saved into a single hdf5_ file. In this module 
hdf5 files are generated and read using h5py_. A typical saved file 
is organised as follows::

	file/
	-- typedict/
	   -- numberOfTypes, (1), uint64
	   -- names, (numberOfTypes), string
	   -- radii, (numberOfTypes), float64
	   -- diffusionConstants, (numberOfTypes), float64
	   -- reactionRadii, (numberOfTypes), float64
	   -- forceTypes, (numberOfTypes), uint64
	-- xyz/
	   -- numberOfParticles, (1), uint64
	   -- typeIds, (numberOfParticles), uint64
	   -- positions, (numberOfParticles, 3), float64
	-- properties/
	   -- temperature, (1), float64
	   -- boxsize, (1), float64
	   -- isReversible, (1), bool
	   -- isPeriodic, (1), bool
	   -- timestep, (1), float64
	   -- repulsionStrength, (1), float64
	-- reactions/ 
		... (to be continued)

``foo/``  marks a group and ``bar, (shape), type`` marks a 
dataset named bar with a given shape and datatype. For more 
information on hdf5 organisation/manipulation, see [h5py_].

.. _hdf5: http://www.hdfgroup.org/HDF5/whatishdf5.html
.. _h5py: http://www.h5py.org

"""
import numpy as np
import sim
import h5py
def saveSimulation(filename, simulation):
	"""
	Save the state of a simulation to file given
	its *simulation* object.

	:param filename: the file will be given this name
	:type filename: string
	:param simulation: object representing your simulation
	:type simulation: pySimulation

	"""
	fileHandle = h5py.File(filename, "w")

	# retrieve typedict info from simulation
	names = simulation.getDictNames()
	radii = simulation.getDictRadii()
	diffs = simulation.getDictDiffusionConstants()
	reactionRadii = simulation.getDictReactionRadii()
	forceTypes = simulation.getDictForceTypes()
	numberOfTypes = len(names)
	# construct typedict group
	typedict = fileHandle.create_group("typedict")
	typedict.create_dataset("numberOfTypes", (1,), dtype=np.uint64)
	typedict["numberOfTypes"][0] = numberOfTypes
	typedict.create_dataset("names", (numberOfTypes,), dtype="S30")
	typedict.create_dataset("radii", (numberOfTypes,), dtype=np.float64)
	typedict.create_dataset(
		"diffusionConstants",
		(numberOfTypes,),
		dtype=np.float64)
	typedict.create_dataset(
		"reactionRadii",
		(numberOfTypes,),
		dtype=np.float64)
	typedict.create_dataset(
		"forceTypes",
		(numberOfTypes,),
		dtype=np.uint64)
	for i in range(numberOfTypes):
		typedict["names"][i]              = names[i]
		typedict["radii"][i]              = radii[i]
		typedict["diffusionConstants"][i] = diffs[i]
		typedict["reactionRadii"][i]      = reactionRadii[i]
		typedict["forceTypes"][i]         = forceTypes[i]

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
	properties["temperature"][0]       = simulation.temperature
	properties["boxsize"][0]           = simulation.boxsize
	properties["isReversible"][0]      = simulation.isReversible
	properties["isPeriodic"][0]        = simulation.isPeriodic
	properties["timestep"][0]          = simulation.timestep
	properties["repulsionStrength"][0] = simulation.repulsionStrength

	# construct reactions group
	# ...
	fileHandle.close()

def loadSimulation(filename):
	"""
	Load the state of simulation from file and return
	an according *simulation* object.

	:param filename: the file to be loaded
	:type filename: string
	:returns: pySimulation object

	"""
	simulation = sim.pySimulation(hasDefaultTypes=False)
	fileHandle = h5py.File(filename, "r")

	# retrieve properties
	properties = fileHandle["properties"]
	simulation.temperature       = properties["temperature"][0]
	simulation.boxsize           = properties["boxsize"][0]
	simulation.isReversible      = properties["isReversible"][0]
	simulation.isPeriodic        = properties["isPeriodic"][0]
	simulation.timestep          = properties["timestep"][0]
	simulation.repulsionStrength = properties["repulsionStrength"][0]

	#retrieve typedict
	typedict = fileHandle["typedict"]
	numberOfTypes = typedict["numberOfTypes"][0]
	for i in range(numberOfTypes):
		simulation.new_Type(
			typedict["names"][i],
			typedict["radii"][i],
			typedict["diffusionConstants"][i],
			typedict["reactionRadii"][i],
			typedict["forceTypes"][i] )

	#retrieve positions
	xyz = fileHandle["xyz"]
	numberOfParticles = xyz["numberOfParticles"][0]
	for i in range(numberOfParticles):
		simulation.addParticle(
			xyz["positions"][i],
			xyz["typeIds"][i] )
	
	#retrieve reactions ...

	return simulation
