##==============================================================================
## montecarlo_driver.jl
##
## Tony Wong (aewsma@rit.edu)
##==============================================================================

cd("/Users/aewsma/codes/CIAM_adaptation_regimes/ciam-code/src")

using Mimi
using MimiCIAM
using RData
using Query
using StatsBase
using NetCDF
using DataFrames

include("montecarlo/ciamMonteCarlo.jl")
include("montecarlo/defmcs.jl")
include("montecarlo/run_ciam_mcs.jl")
include("brickLSL.jl")
include("processResults.jl")


# Define input trial parameters: BRICK model, number of trials (n), min and max BRICK percentile,
#   start and end year, timestep and number of timesteps
trial_params=Dict(
    "brickfile" => "../data/lslr/BRICK_projections.RData",
    "n" => 1000,
    "high" => 97.5,
    "low" => 2.5,
    "ystart" => 2010,
    "yend" => 2150,
    "tstep" => 10,
    "t" => 15,
    "subset" => false,
    "b" => "lsl_rcp85_p50.csv" # won't matter, overwritten when doing the SLR sampling anyway
)

adaptRegime1=Dict(
    "allowMaintain"=>false,
    "noRetreat"=>false,
    "popval"=>1, # 0 = original, 1 = Jones and O'Neill 2016, 2 = Merkens et al 2016
    "fixed"=>true,
    "SSP" => "IIASAGDP_SSP5_v9_130219",
    "SSP_simp" => 5, # won't matter for SSP 0 (old population data) case
    "regimeNum"=>1
)

textheader="run_name,lslr,subset,ssp,ssp_simplified\n"
textstr = "base,$(trial_params["b"]),$(trial_params["subset"]),$(adaptRegime1["SSP"]),$(adaptRegime1["SSP_simp"])"
txtfile=open("../data/batch/init.txt","w") do io
    write(io,textheader)
    write(io,textstr)
end

# vary SLR and CIAM parameters
trial_params["low"] = 2.5
trial_params["high"] = 97.5
runname = "SSP5_BRICK85_varySLR_varyCIAM"
runTrials(85,trial_params,adaptRegime1,vary_slr=true,vary_ciam=true,runname=runname)

# vary only CIAM parameters
trial_params["low"] = trial_params["high"] = 50
runname = "SSP5_BRICK85_varyCIAM"
runTrials(85,trial_params,adaptRegime1,vary_slr=false,vary_ciam=true,runname=runname)

# vary only BRICK parameters
trial_params["low"] = 2.5
trial_params["high"] = 97.5
runname = "SSP5_BRICK85_varyBRICK"
runTrials(85,trial_params,adaptRegime1,vary_slr=true,vary_ciam=false,runname=runname)

# only 5th percentile of SLR, with CIAM defaults
trial_params["low"] = trial_params["high"] = 5
trial_params["n"] = 1
runname = "SSP5_BRICK85_p05"
runTrials(85,trial_params,adaptRegime1,vary_slr=false,vary_ciam=false,runname=runname)

# only 50th percentile of SLR, with CIAM defaults
trial_params["low"] = trial_params["high"] = 50
trial_params["n"] = 1
runname = "SSP5_BRICK85_p50"
runTrials(85,trial_params,adaptRegime1,vary_slr=false,vary_ciam=false,runname=runname)

# only 95th percentile of SLR, with CIAM defaults
trial_params["low"] = trial_params["high"] = 95
trial_params["n"] = 1
runname = "SSP5_BRICK85_p95"
runTrials(85,trial_params,adaptRegime1,vary_slr=false,vary_ciam=false,runname=runname)

##==============================================================================
## End
##==============================================================================
