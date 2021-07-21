# CIAM Uncertainty Propagation

This repository accompanies a study by Wong, Ledna, Rennels, Sheets, Errickson, and Anthoff (2021) using the MimiCIAM implementation of the CIAM model described in [Diaz (2016)](https://doi.org/10.1007/s10584-016-1675-4) ("Original CIAM").

Relative to Original CIAM, MimiCIAM:
* uses updated forcing via the Shared Socioeconomic Pathways (SSPs),
* uses updated sea-level rise projections,
* includes a bug fix that does not allow flood defense heights to be lowered over time,
* and implements a limited decision-maker foresight regarding future sea-level rise.

## Getting started

### Obtain the git repository for this project

First, [clone or download](https://git-scm.com/book/en/v2/Git-Basics-Getting-a-Git-Repository) the `CIAM_uncertainty_propagation` repository. Once this is on your computer, open the Julia console and set this folder as your working directory in the Julia console.

### Installing required packages

This code was created using [Julia v1.5](https://julialang.org/downloads/) and requires several Julia packages.

To install these packages, first enter the package manager by pushing the `]` key in the Julia console. Then, run the following code.

```
add CSV
add DataDeps
add DataFrames
add Dates
add DelimitedFiles
add Distributions
add Mimi
add Missings
add NetCDF
add Query
add RData
add Statistics
add StatsBase
```

You will need to add the `MimiCIAM` model package itself as well. To do this, run the following line while still in the Julia package manager.

```
add https://github.com/BRICK-SLR/MimiCIAM.jl.git
```

Then, press the `Backspace` key to return to Julia.

### Grab the large data files

TODO - down the BRICK projections file

## Run the baseline MimiCIAM simulations

This set of simulations corresponds to Figures 2 and 3 in the manuscript.
* Figure 2 is based on simulations that check the MimiCIAM implementation against the Original CIAM results from the GAMS model version. This includes a version of MimiCIAM that should match the Original CIAM results (correcting for perfect decision-maker foresight for sea-level rise and the construction height decrease fix) and the new "baseline" version of MimiCIAM that uses the Original CIAM forcing data.
* Figure 3 is based on simulations that implement updates to the GDP and population forcing (SSP scenarios), the sea-level rise forcing, and both. Each of these three cases is compared to the baseline MimiCIAM simulation from Figure 2.

To run these baseline comparison simulations, change directories into the `ciam-code/src` directory of the `CIAM_uncertainty_propagation` project. If you were already in the top-level `CIAM_uncertainty_propagation` directory, then you can do this in the Julia console using:
```
cd("ciam-code/src")
```

Then, before running the `baseline_comparisons.jl` script, you will need to change Line 44 where the `brickfile` is set to point to wherever you placed the BRICK sea-level projections locally on your computer. Once you have made that change, run the script either interactively or by running in the Julia console:
```
include("baseline_comparisons.jl")
```

By default, output files will be saved in `ciam-code/output/baseline_comparisons`.

## Run the Monte Carlo ensembles

TODO

## Analysis

TODO - Jupyter notebooks, CSV output

---

Questions? Feedback? Tony Wong (aewsma at rit.edu)
