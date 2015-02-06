from distutils.core import setup
from Cython.Build import cythonize

setup(
	name = "sim",
	ext_modules = cythonize("sim.pyx")
)
