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

# Cite

If you find this library useful please cite us (you might need LaTeX's
`url` package)

    @Misc{Cortes2022,
      author       = {Cortés-Ortuño, David and Fabian, Karl and de Groot, Lennart V.},
      title        = {{micromagnetic_demag_signature}},
      howpublished = {Zenodo \url{doi:10.5281/zenodo.6211355}. Github: \url{https://github.com/Micromagnetic-Tomography/micromagnetic_demag_signature}},
      year         = {2022},
      doi          = {10.5281/zenodo.6211355},
      url          = {https://doi.org/10.5281/zenodo.6211355},
    }

If you have a new version of `biblatex` you can also use `@Software` instead of 
`@Misc`, and add a `version={}` entry.
