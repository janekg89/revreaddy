### TODO
* Fusion reaction that conserves positions of particles that have 0 diffusivity
* trajectory only as .dat so far, make trajectory_unique with .h5,
  therefore extend with hdf5 c-API (only 1 command) and write to extended hyperslab. enable chunking policy
* SimulationImpl child that propagates single moves with a=1
* sample reaction distance from discrete probability distribution. -> new Fusion variant.
* implement clear-period of observables (only timedependant observables) so that they are flushed in predefined intervals and don't fill up memory. Test that buffer is cleared.

#### integration tests:
* periodic boundaries correctly propagated
* minimum image convention fulfilled
* correct force and energy calculation for soft, lj interaction for a system of few particles
* check if reactive partners are found
* check saving/loading oldState
* check correct calculation of acceptance in diffusion and reaction
* check if neighborlist forces yield the same result as naive looping

#### unittests:
* check py wrappers, must be done in python convenience layer
* interactions->calculateForceEnergy and interactions->calculateEnergy
* geometry->doesInteract
