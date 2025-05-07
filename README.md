This repository contains the source files for the dlimIM R package that implements distributed lag interaction models with index modification as described in the paper:

“Distributed Lag Interaction Model with Index Modification” to be submitted to Biostatistics.
Simulation scripts from the manuscript can be found at https://github.com/ddemateis/dlimIM_simulations.

The dlimIM package estimates a distributed lag model with modification by a multiple factors through an index. If you are instead interested in a DLM with modification by a single continuous variable, see the dlim package. If you are instead interested in a DLM with modification by a single categorical or binary variable, see the bdlim package. If you are interested in distributed lag models with heterogeneity by with multiple modifiers, see the heterogeneous distributed lag model in the dlmtree package.

[![DOI](https://zenodo.org/badge/841631784.svg)](https://doi.org/10.5281/zenodo.15354391)

Installation
Install from GitHub:

remotes::install_github("ddemateis/dlimIM")
