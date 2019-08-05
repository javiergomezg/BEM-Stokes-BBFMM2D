from distutils.core import setup
from distutils.extension import Extension

from Cython.Build import cythonize
import os
import numpy
"""
    Creating Module using Cython
"""
setup(
    name="Treecode MPI Wrapper",
    ext_modules = cythonize( #Compile the .pyx to .pyc
        Extension("DL_potential_dir",
                          sources=["src/DL_potential_dir.pyx", "src/H2_2D_Node.cpp", "src/H2_2D_Tree.cpp", "src/kernel_Base.cpp", "src/kernel_Types.cpp"],
                          language = "c++",
                          swig_opts=['-c++','-py3'],
                          include_dirs=[numpy.get_include()],
                          extra_compile_args=['-fPIC', '-O3', '-funroll-loops', '-msse3', '-fopenmp'],
                          extra_link_args=['-fopenmp'],
        )
    )
)
