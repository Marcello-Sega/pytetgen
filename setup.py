# distutils: include_dirs = pytetgen
# distutils: language=c++
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize
import numpy as np
import os 
import codecs

# Get the long description from the README file

here = os.path.abspath(os.path.dirname(__file__))
with codecs.open(os.path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()




setup(	name = 'pytetgen',
	packages=['pytetgen'],
	version = '0.2.1',
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
       		'Development Status :: 4 - Beta',

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
       		'Programming Language :: Python :: 3',
    	],
        ext_modules=[
              Extension('pytetgen', 
                 sources=['pytetgen/tetgen.cxx','pytetgen/predicates.cxx','pytetgen/pytetgen.pyx'],
                 compiler_directives={'language_level':3},
		 include_dirs=[np.get_include()],
		 extra_compile_args=['-g0','-O0','-DTETLIBRARY'],
                 language='c++'),
        ],
        cmdclass = {'build_ext': build_ext},
        install_requires=['triangle>=20200424'],
)
