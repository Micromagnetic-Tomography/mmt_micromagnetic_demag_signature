[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6211355.svg)](https://doi.org/10.5281/zenodo.6211355)

# Micromagnetic Demag Signature

This library allows to calculate the demagnetizing field signature of a
grain onto a two dimensional surface. The grain is modelled micromagnetically
using the MERRILL code for the simulation of magnetic samples using finite
element discretization.

## Theory

The demagnetizing field signal is modelled by taking every magnetic moment from
the nodes of the finite element mesh and calculating its dipolar signal in
the two dimensional surface.
