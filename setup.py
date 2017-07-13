from distutils.core import setup,Extension
from Cython.Build import cythonize
import numpy as np
setup(	name = 'pytetgen',
	version = '0.1',
	description = 'wrapper for the tetgen mesh generator',
	author = 'Marcello Sega',
	author_email = 'marcello.sega@gmail.com',
	ext_modules = cythonize(
                       "pytetgen.pyx",                 # our Cython source
                       sources=["tetgen.cxx","predicates.cxx"],  # additional source file(s)
                       language="c++",             # generate C++ code
                    ),
     	include_dirs = [np.get_include()]
)
