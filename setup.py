from distutils.core import setup, Extension

module1 = Extension('revreaddy',
			include_dirs = ['/usr/local/include',
					'/usr/include',
					'include/'],
			libraries=['gsl', 'gslcblas', 'm'], 
			library_dirs = ['/usr/local/lib'],
			sources = 	['src/revreaddymodule.cpp',
						'src/Random.cpp',
						'src/Particle.cpp',
						'src/ActiveParticles.cpp',
						'src/Simulation.cpp',
						'src/SingleParticleDiffusion.cpp'],
			extra_compile_args = ['-std=c++11'])

setup(name='revreaddy', 
		version='0.0',
		description='reversible reaction diffusion dynamics simulation',
		ext_modules=[module1])
