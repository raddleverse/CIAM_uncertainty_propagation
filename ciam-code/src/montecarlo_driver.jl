using DataFrames: Highlighter
##==============================================================================
## montecarlo_driver.jl
##
## Tony Wong (aewsma@rit.edu)
##==============================================================================

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

brickfile = "/Users/lisarennels/JuliaProjects/CIAMPaper/local-data/BRICK_projections.RData"
outputdir = joinpath(@__DIR__, "..", "output", "MonteCarlo")
isdir(outputdir) || mkdir(outputdir)

# write the init file
init_settings = Dict(
    :init_filename   => "MCdriver_init.csv",
    :subset          => false,
    :ssp             => "IIASAGDP_SSP5_v9_130219",
    :ssp_simplified  => 5 # won't matter, just to get defaults for the popdens_seg_jones and _merkens arrays.
    :b => "lsl_rcp85_p50.csv" # won't matter, overwritten when doing the SLR sampling anyway
)

init_file = joinpath(outputdir,init_settings[:init_filename])
textheader="run_name,lslr,subset,ssp,ssp_simplified\n"
textstr = "base,$(init_settings[:b]), $(init_settings[:subset]), $(init_settings[:ssp]), $(init_settings["ssp_simplified"])"
txtfile=open(init_file,"w") do io
    write(io,textheader)
    write(io,textstr)
end

# Define input trial parameters: BRICK model, number of trials (n), min and max BRICK percentile,
# start and end year, and timestep

trial_params = Dict(
    :brickfile  => brickfile,
    :n      => 1000,
    :high   => 97.5,
    :low    =>  2.5,
    :ystart => 2010,
    :yend   => 2150,
    :tstep  => 10,
)

# Define other eneded parameters and settings for the model
adaptRegime1 = Dict(
    :fixed          => true,
    :t              => 15,
    :noRetreat      => false,
    :allowMaintain  => false,
    :popval         => 1,
    :GAMSmatch      => true,
    :subset         => false
)

# vary SLR and CIAM parameters
trial_params[:low] = 2.5
trial_params[:high] = 97.5
runname = "SSP5_BRICK85_varySLR_varyCIAM"
runTrials(85, trial_params, adaptRegime1, outputdir, init_file, vary_slr=true, vary_ciam=true, runname=runname)

# vary only CIAM parameters
trial_params[:low] = 50
trial_params[:high] = 50
runname = "SSP5_BRICK85_varyCIAM"
runTrials(85, trial_params, adaptRegime1, outputdir, init_file, vary_slr=false, vary_ciam=true, runname=runname)

# vary only BRICK parameters
trial_params[:low] = 2.5
trial_params[:high] = 97.5
runname = "SSP5_BRICK85_varyBRICK"
runTrials(85, trial_params, adaptRegime1, outputdir, init_file, vary_slr=true, vary_ciam=false, runname=runname)

# only 5th percentile of SLR, with CIAM defaults
trial_params[:low] = 5
trial_params[:high] = 5
trial_params[:n] = 1
runname = "SSP5_BRICK85_p05"
runTrials(85, trial_params, adaptRegime1, outputdir, init_file, vary_slr=false, vary_ciam=false, runname=runname)

# only 50th percentile of SLR, with CIAM defaults
trial_params[:low] = 50
trial_params[:high] = 50
trial_params[:n] = 1
runname = "SSP5_BRICK85_p50"
runTrials(85, trial_params, adaptRegime1, outputdir, init_file, vary_slr=false, vary_ciam=false, runname=runname)

# only 95th percentile of SLR, with CIAM defaults
trial_params[:low] = 95
trial_params[:high] = 95
trial_params[:n] = 1
runname = "SSP5_BRICK85_p95"
runTrials(85, trial_params, adaptRegime1, outputdir, init_file, vary_slr=false, vary_ciam=false, runname=runname)

##==============================================================================
## End
##==============================================================================
