from distutils.core import setup
from Cython.Build import cythonize
import numpy

setup(
    ext_modules=cythonize(["satfit_accelerated.pyx"]),
	include_dirs=[numpy.get_include()],
    extra_compile_args=["-O3","-ffast-math"]
)

