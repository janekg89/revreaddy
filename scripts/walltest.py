import time
import numpy as np
import revreaddy.sim as sim

simulation = sim.pySimulation()
simulation.boxsize = 20.
simulation.isPeriodic = True
simulation.isReversible = False
simulation.temperature = 1.
simulation.repulsionStrength = 1.

simulation.new_Type("wall-particle", 1.0, 1.0, 0.0, 2)
simulation.new_Wall([0.,0.,-1.], [0.,0., 1.], 1., [3])
simulation.new_Wall([0.,0., 1.], [0.,0.,-1.], 1., [3])

L = simulation.boxsize / 2.

for x in np.arange(-L,L,2.):
	for y in np.arange(-L,L,2.):
		position = [x, y, 0.]
		simulation.addParticle(position, 3)

trajname = "data/traj_" + time.strftime("%Y_%m_%d-%H_%M_%S") + ".xyz"
simulation.new_Trajectory(trajname)
simulation.maxTime = 20000
simulation.timestep = 1e-4
simulation.run()

simulation.writeAllObservablesToFile()
simulation.deleteAllObservables()
