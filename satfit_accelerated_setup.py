from distutils.core import setup
from Cython.Build import cythonize
from Cython.Distutils import build_ext
import numpy

setup(
    ext_modules=cythonize(["satfit_caccelerated.pyx"]),
	include_dirs=[numpy.get_include()],
    extra_compile_args=["-O3","-ffast-math"]
)

