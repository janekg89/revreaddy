revreaddy v0.1
**************

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
cython. Large parts are still undocumented. 
Not many unittests.

Installation
============

Requirements for compiling the cython extension:
	* GNU scientific libraries (libgsl)
	  for random number generation
	* HDF5 (libhdf5/libhdf5_hl) for saving observables 
	  to binary file formats

Requirements for usage:
	* Python 2
	* h5py python module is required for saving
	  and loading simulations with the
	  utils.py module

Build the extension in place with

::

	$ python setupCython.py build_ext --inplace

This will generate the sim.so file, which is basically
the whole simulation code. When building in place it
is convenient to add the parent directory of revreaddy
to your PYTHONPATH by typing

::

	$ export PYTHONPATH=/parentPathOfRevreaddy/

in a bash terminal or add this line to your .bashrc file.
You can now start working by importing revreaddy in a
python script, e.g.

::

	>>> import revreaddy.sim as sim
	>>> import revreaddy.utils as utils

Unittests
=========

As a unittest framework, cxxtest is used on the c-side.
These tests are build with

::

	$ make unittest

and run with

::

	$ bin/unittest

