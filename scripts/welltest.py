import time
import numpy as np
import revreaddy.sim as sim

simulation = sim.pySimulation()
simulation.boxsize = 30.
simulation.isPeriodic = True
simulation.isReversible = False
simulation.temperature = 1.
simulation.repulsionStrength = 0.

simulation.new_Type("well-particle", 1.0, 1.0, 0.0, 2)
distanceMinima = 0.15 * simulation.boxsize
wellStrength = 1.
simulation.new_DoubleWellZ(distanceMinima, wellStrength, [3])

L = simulation.boxsize / 2.

#for x in np.arange(-L,L,2.):
#	for y in np.arange(-L,L,2.):
#		position = [x, y, 0.]
position = [0., 0., 0.]
simulation.addParticle(position, 3)

simulation.maxTime = 2000000
simulation.run()

#trajname = "data/traj_" + time.strftime("%Y_%m_%d-%H_%M_%S") + ".xyz"
probname = "data/prob_" + time.strftime("%Y_%m_%d-%H_%M_%S") + ".dat"
#simulation.new_Trajectory(trajname)
ranges = np.arange(-10., 10., 0.2)
simulation.new_ProbabilityDensity(probname, 3, ranges, 2)
simulation.maxTime = 20000000
simulation.timestep = 1e-3
simulation.run()

simulation.writeAllObservablesToFile()
simulation.deleteAllObservables()

import matplotlib.pyplot as plt
prob = np.loadtxt(probname)
plt.plot(prob[0], prob[1])
plt.show()
