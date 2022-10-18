[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6211355.svg)](https://doi.org/10.5281/zenodo.6211355)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

# MMT Numerical Libraries: Micromagnetic Demag Signature

This library allows to calculate the demagnetizing field signature of a
grain onto a two dimensional surface. The grain is modelled micromagnetically
using the MERRILL code for the simulation of magnetic samples using finite
element discretization.

This library is developed as part of the ![Mimatom / MMT](https://mimatom.org/)
project.

## Theory

The demagnetizing field signal is modelled by taking every magnetic moment from
the nodes of the finite element mesh and calculating its dipolar signal in the
two dimensional surface. 

## Library

The magnetic moments are computed from the magnetization field and node volumes
of the mesh nodes which are obtained from MERRILL's `vbox` output files. These
files can be blank space separated (old MERRILL versions) or comma separated.
The magnetic field calculations are parallelized using `cython`.

A tutorial is planned for future releases. Refer to the `test` folder to find
examples using this library.

# Cite

If you find this library useful please cite us (you might need LaTeX's
`url` package)

    @Misc{Cortes2022,
      author       = {Cortés-Ortuño, David and Fabian, Karl and de Groot, Lennart V.},
      title        = {{MMT Numerical Libraries: Micromagnetic Demag Signature}},
      publisher    = {Zenodo},
      note         = {Github: \url{https://github.com/Micromagnetic-Tomography/mmt_micromagnetic_demag_signature}},
      year         = {2022},
      doi          = {10.5281/zenodo.6211355},
      url          = {https://doi.org/10.5281/zenodo.6211355},
    }

If you have a new version of `biblatex` you can also use `@Software` instead of 
`@Misc`, and add a `version={}` entry.
