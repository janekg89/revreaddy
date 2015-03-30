import numpy as np
import matplotlib.pyplot as plt
import time
import sim
simulation = sim.pySimulation()
simulation.isPeriodic        = True
simulation.boxsize           = 10.
simulation.timestep          = 0.0001
simulation.temperature       = 1.
simulation.repulsionStrength = 1.
simulation.verbose           = False

particleType      = simulation.particleTypeToId("soft")
radius            = 1.
diffusionConstant = 1.

L = simulation.boxsize / 2.
for x in np.arange(-L, L, 2.*radius):
    for y in np.arange(-L, L, 2.*radius):
        for z in np.arange(-L, L, 2.*radius):
            position = np.array([x, y, z])
            simulation.addParticle(
                position,
                particleType,
                radius,
                diffusionConstant
            )

simulation.new_Trajectory("traj.xyz")

t1 = time.clock()
simulation.maxTime = 10
simulation.run()
t2 = time.clock()

simulation.writeAllObservablesToFile()

print "Time for", simulation.getParticleNumber(),"particles and", \
      simulation.maxTime ,"timesteps:", t2 - t1, "seconds"