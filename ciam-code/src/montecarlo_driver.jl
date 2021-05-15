# todo - monte carlo simulations with
# 1) updated SLR projections (BRICK)

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
    "n" => 10,
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
    "popval"=>1,
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

runTrials(85,trial_params,adaptRegime1)
