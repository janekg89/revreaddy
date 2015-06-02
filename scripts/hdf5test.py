import revreaddy.sim as sim
import revreaddy.utils as utils

simulation = sim.pySimulation()
simulation.boxsize = 12.
simulation.timestep = 1e-5
simulation.temperature = 2.
simulation.addParticle([0.,0.,0.], 2)
simulation.addParticle([0.,0.,1.], 2)
simulation.addParticle([0.,1.,0.], 2)
simulation.addParticle([1.,0.,0.], 2)
simulation.addParticle([0.,0.,2.], 1)
simulation.addParticle([0.,2.,0.], 1)

simulation.addParticle([-1.,-1.,-1.], 3)
simulation.new_Type("H2O", 2.3, 3.14, 5., 1)
simulation.addParticle([-1.,-1.,-1.], 3)
#utils.saveSimulation("ubertest.hdf5", simulation)
