import setuptools
from setuptools.extension import Extension
# import sys
# cython and python dependency is handled by pyproject.toml
from Cython.Build import cythonize
import numpy


# -----------------------------------------------------------------------------
# Compilation of C module in c_lib
com_args = ['-std=c99', '-O3', '-fopenmp']
link_args = ['-fopenmp']
extensions = [
    Extension("micromagnetic_demag_signature.clib.mds_clib",
              ["micromagnetic_demag_signature/clib/mds_clib.pyx",
               "micromagnetic_demag_signature/clib/clib.c"],
              extra_compile_args=com_args,
              extra_link_args=link_args,
              include_dirs=[numpy.get_include()]
              ),
]

# -----------------------------------------------------------------------------

with open('README.md') as f:
    long_description = f.read()

setuptools.setup(
    # setup_requires=['cython'],  # not working (see the link at top)
    name='micromagnetic_demag_signature',
    version='0.1',
    description=('Python lib to calculate the signature of the demag field '
                 'from a micromagnetic model of a grain, onto a 2D surface'),
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='D. Cortés-Ortuño, K. Fabian, L. V. de Groot',
    author_email='d.i.cortes@uu.nl',
    packages=setuptools.find_packages(),
    ext_modules=cythonize(extensions),
    install_requires=['matplotlib',
                      'numpy>=1.20',
                      'scipy>=1.6',
                      'numba>=0.51',
                      ],

    # TODO: Update license
    classifiers=['License :: BSD2 License',
                 'Programming Language :: Python :: 3 :: Only',
                 ],
)
