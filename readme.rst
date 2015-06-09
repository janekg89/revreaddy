"""

A particle based reaction-diffusion simulation with a
reversible integrator, which obeys detailed balance.

This early version features:
	* Brownian Dynamics integrator with Metropolis-Hastings
	  correction
	* NO reactions (not yet!)
	* customizable particle types with different radii and
	  diffusion constants
	* two possible repulsion potentials, harmonic repulsion
	  and Lennard-Jones interaction
	* geometric building blocks implemented as single
	  particle potentials, currently repulsive plane "wall" 
	* cubic periodic boundary conditions, can be switched
	  off (not recommended) 
	* observables which are calculated on the fly and saved
	  to binary (.hdf5/.h5) files or text files
	  (.dat/.txt). These are acceptance, energy,
	  mean-squared-displacement, probability density,
	  radial distribution function, trajectory
	* the current state of the simulation can be saved
	  to a compact binary format with a single line
	  and also loaded from file with a single line

revreaddy is mainly a python module which wraps the
functionality of an underlying C++ software using
cython.

Large parts are still undocumented and unittests are not
up-to-date.

"""
