## revreaddy

A particle based reaction-diffusion simulation with a
reversible integrator, which obeys detailed balance.
Revreaddy is mainly a C++ software which is wrapped to
be accessible from python. Main purpose of this
software is to try out different reaction-diffusion
schemes.

This version features:
* Brownian Dynamics integrator with optional Metropolis-Hastings
  correction
* Reactions that allow particle creation/destruction
  e.g. Fusion reaction A + B <--> C,
  with forward and backward rates
* customizable particle types with different radii and
  diffusion constants
* two possible repulsion potentials, harmonic repulsion
  and Lennard-Jones interaction
* geometric building blocks implemented as single
  particle potentials, currently repulsive plane "wall"
  and a double-well potential in one dimension
* cubic periodic boundary conditions, can be switched
  off (not recommended) 
* observables which are calculated on the fly and saved
  to binary (.hdf5/.h5) files or text files
  (.dat/.txt). These are acceptance, energy,
  mean-squared-displacement, probability density,
  radial distribution function, trajectory
* Neighborlattice force calculation leading to overall
  complexity *O(#particles)* (exception when rdf is
  calculated on the fly)

### Dependencies

Linked libraries:
* GNU scientific libraries (gsl, gslcblas)
  for random number generation and histogram
* HDF5 (hdf5, hdf5_hl, H5Cpp) for saving observables 
  to binary file formats
* Boost for python wrapping, logging and stuff
  (Boost.python, Boost.log, Boost.system, Boost.filesystem,
  Boost.multi_array)

Requirements for usage:
* Python
* h5py

### Building
The way to build this using cmake goes like

    $ mkdir build
    $ cd build
    $ cmake .. $CMAKE_FLAGS
    $ make -j
    $ cd ..

Any user specific configuration goes into the `$CMAKE_FLAGS` 
and should be defined beforehand.

### Usage
Import revreaddy in a python script, e.g.

    >>> import revreaddy 
    >>> s = revreaddy.Sim()
    >>> s.run(steps=10000, timestep=0.1)