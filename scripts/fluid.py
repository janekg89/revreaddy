import numpy as np
import time
import sim
simulation = sim.pySimulation()
simulation.isPeriodic        = True
simulation.boxsize           = 10.
simulation.timestep          = 0.0001
simulation.temperature       = 1.
simulation.repulsionStrength = 1.
simulation.verbose           = False

particleType      = sim.typeStringToId['soft']

L = simulation.boxsize / 2.
for x in np.arange(-L, L, 2.):
    for y in np.arange(-L, L, 2.):
        for z in np.arange(-L, L, 2.):
            position = np.array([x, y, z])
            simulation.addParticle(
                position,
                particleType
            )

def go():
	t1 = time.clock()
	simulation.maxTime = 3000
	simulation.run()
	t2 = time.clock()

	print "Time for", simulation.getParticleNumber(),"particles and", \
      simulation.maxTime ,"timesteps:", t2 - t1, "seconds"
