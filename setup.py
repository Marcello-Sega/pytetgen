# distutils: include_dirs = tetgen
from distutils.core import setup,Extension
from Cython.Build import cythonize
import numpy as np
import os 
import codecs

# Get the long description from the README file

here = os.path.abspath(os.path.dirname(__file__))
with codecs.open(os.path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

setup(	name = 'pytetgen',
	version = '0.1.0',
	description = 'wrapper for the tetgen mesh generator',
	long_description=long_description,
	author = 'Marcello Sega',
	author_email = 'marcello.sega@gmail.com',
        url='https://github.com/Marcello-Sega/pytetgen',
   	license='AGPLv3',
	classifiers=[
       		# How mature is this project? Common values are
       		#   3 - Alpha
       		#   4 - Beta
       		#   5 - Production/Stable
       		'Development Status :: 3 - Alpha',

       		# Indicate who your project is intended for
       		'Intended Audience :: Science/Research',
       		'Topic :: Scientific/Engineering :: Bio-Informatics',
       		'Topic :: Scientific/Engineering :: Chemistry',
       		'Topic :: Scientific/Engineering :: Physics',
       		'Topic :: Software Development :: Libraries :: Python Modules',

       		# Pick your license as you wish (should match "license" above)
   		'License :: OSI Approved :: GNU Affero General Public License v3 or later (AGPLv3+)',

       		# Specify the Python versions you support here. In particular, ensure
       		# that you indicate whether you support Python 2, Python 3 or both.
       		'Programming Language :: Python :: 2.7',
    	],

	ext_modules = cythonize(
                       "pytetgen/pytetgen.pyx",                 # our Cython source
                       sources=["tetgen/tetgen.cxx","tetgen/predicates.cxx"],  # additional source file(s)
                       language="c++",             # generate C++ code
                    ),
     	include_dirs = ['./',np.get_include(),'tetgen']
)