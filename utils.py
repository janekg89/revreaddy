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
	-- forcemap/
	   -- numberOfForces, (1), uint64
	   -- names, (numberOfForces), string
	   -- types, (numberOfForces), string
	   -- affectedTuples, (numberOfForces, 2), uint64
	   -- parameters, (numberOfForces, 3), float64
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
	#forceTypes = simulation.getDictForceTypes()
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
	#typedict.create_dataset(
	#	"forceTypes",
	#	(numberOfTypes,),
	#	dtype=np.uint64)
	for i in range(numberOfTypes):
		typedict["names"][i]              = names[i]
		typedict["radii"][i]              = radii[i]
		typedict["diffusionConstants"][i] = diffs[i]
		typedict["reactionRadii"][i]      = reactionRadii[i]
		#typedict["forceTypes"][i]         = forceTypes[i]

	# construct forcemap group
	forcemap = fileHandle.create_group("forcemap")
	numberOfForces = simulation.getNumberForces()
	forcemap.create_dataset("numberOfForces", (1,), dtype=np.uint64)
	forcemap.create_dataset("names", (numberOfForces,), dtype="S30")
	forcemap.create_dataset("types", (numberOfForces,), dtype="S30")
	forcemap.create_dataset("affectedTuples", (numberOfForces,2), dtype=np.uint64)
	forcemap.create_dataset("parameters", (numberOfForces,3), dtype=np.float64)
	forcemap["numberOfForces"][0] = numberOfForces
	# fill in datasets
	for i in range(numberOfForces):
		forcemap["names"][i] = simulation.getForceName(i)
		forcemap["types"][i] = simulation.getForceType(i)
		forcemap["affectedTuples"][i,0] = simulation.getForceAffectedTuple(i)[0]
		forcemap["affectedTuples"][i,1] = simulation.getForceAffectedTuple(i)[1]
		# if softrepulsion, the only parameter is repulsion strength
		if (simulation.getForceType(i) == "SoftRepulsion"):
			forcemap["parameters"][i,0] = simulation.getForceParameters(i)[0]
		# if lennardjones, the only parameter is epsilon
		elif (simulation.getForceType(i) == "LennardJones"):
			forcemap["parameters"][i,0] = simulation.getForceParameters(i)[0]
		else:
			print "Error: Unknown force type in saveSimulation. Saving aborted"
			fileHandle.close()
			return

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
	#properties.create_dataset("repulsionStrength", (1,), dtype=np.float64)
	properties["temperature"][0]       = simulation.temperature
	properties["boxsize"][0]           = simulation.boxsize
	properties["isReversible"][0]      = simulation.isReversible
	properties["isPeriodic"][0]        = simulation.isPeriodic
	properties["timestep"][0]          = simulation.timestep
	#properties["repulsionStrength"][0] = simulation.repulsionStrength

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
	simulation = sim.pySimulation()#hasDefaultTypes=False)
	fileHandle = h5py.File(filename, "r")

	# retrieve properties
	properties = fileHandle["properties"]
	simulation.temperature       = properties["temperature"][0]
	simulation.boxsize           = properties["boxsize"][0]
	simulation.isReversible      = properties["isReversible"][0]
	simulation.isPeriodic        = properties["isPeriodic"][0]
	simulation.timestep          = properties["timestep"][0]
	#simulation.repulsionStrength = properties["repulsionStrength"][0]

	#retrieve typedict
	typedict = fileHandle["typedict"]
	numberOfTypes = typedict["numberOfTypes"][0]
	for i in range(numberOfTypes):
		simulation.new_Type(
			typedict["names"][i],
			typedict["radii"][i],
			typedict["diffusionConstants"][i],
			typedict["reactionRadii"][i])
			#typedict["forceTypes"][i] )

	#retrieve positions
	xyz = fileHandle["xyz"]
	numberOfParticles = xyz["numberOfParticles"][0]
	for i in range(numberOfParticles):
		simulation.addParticle(
			xyz["positions"][i],
			xyz["typeIds"][i] )

	#retrieve forcemap
	forcemap = fileHandle["forcemap"]
	numberOfForces = forcemap["numberOfForces"][0]
	for i in range(numberOfForces):
		if (forcemap["types"][i] == "SoftRepulsion"):
			simulation.new_SoftRepulsion(
				forcemap["names"][i],
				[forcemap["affectedTuples"][i,0],forcemap["affectedTuples"][i,1]],
				forcemap["parameters"][i,0])

		elif (forcemap["types"][i] == "LennardJones"):
			simulation.new_LennardJones(
				forcemap["names"][i],
				[forcemap["affectedTuples"][i,0],forcemap["affectedTuples"][i,1]],
				forcemap["parameters"][i,0])
		else:
			print "Error: Unknown force type in loadSimulation. Loading aborted. Returning empty simulation"
			fileHandle.close()
			return sim.pySimulation()

	#retrieve reactions ...

	return simulation

def fillCuboidWithParticles(simulation, r1, r2, numberOfParticles, particleType):
	x1, x2 = r1[0], r2[0]
	y1, y2 = r1[1], r2[1]
	z1, z2 = r1[2], r2[2]
	distanceX = (x2 - x1) / float(numberOfParticles)**(1./3.)
	distanceY = (y2 - y1) / float(numberOfParticles)**(1./3.)
	distanceZ = (z2 - z1) / float(numberOfParticles)**(1./3.)
	xRange = np.arange(x1, x2, distanceX) + distanceX / 2.
	yRange = np.arange(y1, y2, distanceY) + distanceY / 2.
	zRange = np.arange(z1, z2, distanceZ) + distanceZ / 2.
	counter = 0
	for x in xRange:
		for y in yRange:
			for z in zRange:
				if (counter < numberOfParticles):
					simulation.addParticle([x,y,z], particleType)
	return

def showSnapshotMatPlotlib(simulation):
	import matplotlib.pyplot as plt
	from mpl_toolkits.mplot3d import Axes3D
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	for i in range(simulation.getParticleNumber()):
		position = simulation.getPosition(i)
		ax.scatter(
			position[0], position[1], position[2],
			color="blue", marker="o", alpha=0.3, s=50.)

def generateSnapshotVmd(filename, simulation):
	N = simulation.getParticleNumber()
	xyzFile = open(filename, "w")
	xyzFile.write(str(N) + "\n")
	xyzFile.write("#\n")
	for i in range(N):
		particleType = str( simulation.getTypeId(i) )
		position = simulation.getPosition(i)
		x = str( position[0] )
		y = str( position[1] ) 
		z = str( position[2] )
		line = "T" + particleType + "\t" + x + "\t" + y + "\t" + z + "\n"
		xyzFile.write(line)
	xyzFile.close()

	# mostly adapted from Johannes Schoeneberg's ReaDDy software
	# see github.com/readdy
	tclScript = open(filename + ".tcl", "w")
	tclScript.write("mol delete top\n")
	tclScript.write("mol load xyz " + filename + "\n")
	tclScript.write("mol delrep 0 top\n")
	tclScript.write("display resetview\n")
	dictRadii = simulation.getDictRadii()
	M = len( dictRadii )
	for i in range(M):
		tclScript.write(
			"mol representation VDW " + str(dictRadii[i] * 0.7) + " 16.0\n"
		)
		tclScript.write("mol selection name T" + str(i) + "\n")
		tclScript.write("mol addrep top\n")
	tclScript.write("animate goto 0\n")
	tclScript.write("color Display Background white\n")
	tclScript.write(
	"molinfo top set {center_matrix} {{{1 0 0 0}{0 1 0 0}{0 0 1 0}{0 0 0 1}}}\n"
	)
	tclScript.close()
