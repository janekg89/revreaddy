from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension
import os

# locate and collect all sourcefiles
src = ["sim.pyx"]
for f in os.listdir( os.path.join(os.getcwd(), "src") ):
	if ( f.endswith(".cpp") and ( f != "main.cpp") ):
		src += [ os.path.join("src", f) ]

print "sources:", src

extensions = [
	Extension(
		"sim",
		sources = src,
		language = "c++",
		extra_compile_args = ["-std=c++11","-O3"],
		include_dirs = ["/usr/include","/usr/local/include","include"],
		libraries = ["m","gsl","gslcblas","hdf5","hdf5_hl"]
	)
]

setup(
	name = "sim",
	ext_modules = cythonize(
		extensions,
		compiler_directives={"embedsignature" : True}
	)
)
