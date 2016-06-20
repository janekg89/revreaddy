### TODO
* Fusion reaction that conserves positions of particles that have 0 diffusivity
* SimulationImpl child that propagates single moves with a=1
* revreaddy.py ->  show_simulation().

#### integration tests:
* correct force and energy calculation for soft, lj interaction for a system of few particles
* check if reactive partners are found
* check saving/loading oldState
* check correct calculation of acceptance in diffusion and reaction

#### unittests:
* check py wrappers, must be done in python convenience layer
* interactions->calculateForceEnergy and interactions->calculateEnergy
* geometry->doesInteract
