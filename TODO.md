### TODO
* implement clear-period of observables so that they are flushed in predefined intervals and don't fill up memory

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