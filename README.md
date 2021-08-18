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
add https://github.com/SLURMverse/MimiCIAM.jl.git
```

Then, press the `Backspace` key to return to Julia.

### Grab the large data files

#### Sea-level projections

If you obtained this repository via Zenodo, then the large data file containing BRICK sea-level projections from [Vega-Westhoff et al. (2019)](https://doi.org/10.1029/2018EF001082) is in this repository at `ciam-code/data/lslr/BRICK_projections.RData`.

If you obtained this repository via GitHub, then this large data file was not present. You can download it by running the following code from the Julia console, in the top-level directory of this project.

```
if !isfile("./ciam-code/data/lslr/BRICK_projections.RData")
    url = "https://zenodo.org/record/3628215/files/sample_projections.RData"
    download(url,  "./ciam-code/data/lslr/BRICK_projections.RData")
end
```

#### Original CIAM GAMS outputs for benchmarking

You will need a directory full of results from the original implementation of CIAM in GAMS from [Diaz (2016)](https://doi.org/10.1007/s10584-016-1675-4) in order to complete the baseline comparisons (Figure 2 in the manuscript). If you obtained this repository via Zenodo, then those files should be in `ciam-code/output/originalCIAM_gams_outputs`.

If you obtained this repository via GitHub, then these instructions will later be modified to link to those files on Zenodo.

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

By default, output files will be saved in `ciam-code/output/baseline_comparisons`. Analyzing this output and generating figures/output data files is described below under **Analysis**.

## Run the Monte Carlo ensembles

To run the Monte Carlo ensembles under SSP5-RCP8.5 (main text) and SSP1-RCP2.6 (Supporting Material), change directories into the `ciam-code/src` directory of the `CIAM_uncertainty_propagation` project. Then, run the `montecarlo_driver.jl` script. It will give you an error about missing the `BRICK_projections.RData` file if you do not have it in the `ciam-code/data/lslr` directory.

The script is set up to run by default ensembles of 1000 simulations, varying sea-level rise scenarios, varying CIAM socioeconomic parameters, and varying both. This is done for each of SSP5-RCP8.5 and SSP1-RCP2.6.

The script will also run single simulations using the sea-level rise ensemble members that give the 5th, 50th, and 95th percentiles of global mean sea-level rise in 2150. CIAM parameters are held fixed at their default central values. The point here is to examine how the anticipated coastal adaptation costs differ between:
* Case 1: you compute these percentiles of the sea-level rise ensemble, and assume they are representative of the same percentiles of the distribution of adaptation costs; and
* Case 2: you compute the adaptation costs for _all_ sea-level rise ensemble members, then compute the percentiles of the distribution of adaptation costs.

## Analysis

To run the analysis for the baseline comparisons, run through the Jupyter notebook that is at `work_baseline_comparisons/plotsAndAnalysis_baselineComparisons.ipynb`. As the notebook mentions, some of those validation steps where the MimiCIAM results are compared against the Original CIAM GAMS results will take a long time. For example, changing the MimiCIAM segments to adjust for perfect foresight (so they'll match the old GAMS results) takes about 30 hours on a desktop workstation. (In the MimiCIAM code moving forward, only a random subset of segments are checked, so it is much faster.)

To run the analysis for the Monte Carlo uncertainty propagation experiments, run through the Jupyter notebooks that are at `work_uncertainty_propagation/SSP5-RCP85.ipynb` and `work_uncertainty_propagation/SSP1-RCP26.ipynb`. These are much faster and should not take more than a few minutes to run interactively.

Both Jupyter notebooks will yield CSV files containing the numbers used to generate figures.

---

Questions? Feedback? Tony Wong (aewsma at rit.edu) We are also interested in helping folks use other sea-level projections with MimiCIAM. So, please also feel free to ask questions and/or open issues in the [MimiCIAM.jl](https://github.com/SLURMverse/MimiCIAM.jl) GitHub repository.
