from setuptools import Extension, setup, find_packages
from Cython.Build import cythonize
import numpy as np

# https://setuptools.readthedocs.io/en/latest/setuptools.html#declaring-extras-optional-features-with-their-own-dependencies

extensions = [
  Extension('trusat.caccelerated',
            sources=['trusat/caccelerated.pyx'],
            include_dirs=[np.get_include()],
            extra_compile_args=['-O3','-ffast-math','-march=native'],
  ),
  Extension('trusat.profile',
            sources=['trusat/profile.pyx'],           
            optional='optional',
            include_dirs=[np.get_include(), '/usr/local/include'],
            libraries=["m"],
            extra_compile_args=['-O3','-ffast-math','-march=native','-fopenmp','-I/usr/local/opt/llvm/include'],
            extra_link_args=['-lomp',"-L/usr/local/opt/libomp/lib/", "-L/usr/local/opt/llvm/lib"],
            # extra_link_args=['-lomp','-L/usr/local/opt/libomp/lib/'],
  )
]

with open("README.md", "r") as fh:
    long_description = fh.read()

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setup(
    name="trusat-orbit", 
    version="0.9.0",
    author="Chris Lewicki",
    author_email="chris@lewicki.com",
    description="TruSat satellite observation processing utilities",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://trusat.org/",
    project_urls={
        "Git repo" : "https://github.com/TruSat/trusat-orbit",
        "Learning Hub" : "https://learn.trusat.org/docs/start-here",
        "Forums" : "https://discuss.trusat.org/"
    },
    packages=['trusat'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=[
        requirements,
        'trusat_backend==1.1.0',
    ],
    ext_modules=cythonize(extensions),
    dependency_links = [
        'git+https://github.com/TruSat/trusat-backend@dev.chris.package#egg=trusat_backend-1.1.0'
    ],
)