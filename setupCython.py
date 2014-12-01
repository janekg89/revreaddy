from distutils.core import setup
from Cython.Build import cythonize

setup(ext_modules = cythonize(
	"sim.pyx",
	sources=["src/Simulation.cpp", "src/Particle.cpp"],
	include_dirs = ['/usr/local/include',
		'/usr/include',
		'include/'],
	libraries=['m'],
	language="C++"
))
