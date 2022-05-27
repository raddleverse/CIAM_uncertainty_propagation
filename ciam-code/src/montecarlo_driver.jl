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

# The file structure created by this process will look as follows:
#
# - top directory: output/MonteCarlo holding subdirectories and the init file MCdriver_init.csv
#
# - for each call to runTrials, a (unique) subdir created with the code:
#   joinpath(outputdir, runname, "CIAM $(Dates.format(now(), "yyyy-mm-dd HH-MM-SS")) MC$(trial_params[:n])")
#
# - for each of these subdirs, you will see a RawResults folder with any results directly
#   from the monte carlo runs in run_ciam_mcs.jl and a PostProcessing folder written
#   with the code at the bottom of ciamMonteCarlo.jl
#

# First, set some things up:

#brickfile = joinpath(@__DIR__, "..", "data", "lslr", "BRICK_projections.RData")
brickfile = "https://zenodo.org/record/6461560/files/sneasybrick_projections_csv.zip"
#brickfile = joinpath(@__DIR__, "..", "data", "lslr", "sneasybrick_projections_csv.zip")
outputdir = joinpath(@__DIR__, "..", "output", "MonteCarlo")
isdir(outputdir) || mkpath(outputdir)

ssp_files = Dict(1 => "IIASAGDP_SSP1_v9_130219",
                 2 => "IIASAGDP_SSP2_v9_130219",
                 3 => "IIASAGDP_SSP3_v9_130219",
                 4 => "IIASAGDP_SSP4_v9_130219",
                 5 => "IIASAGDP_SSP5_v9_130219")
popinput = 0                          # population density input data (only 0 is supported currently)
ssp_rcp_scenarios = [(1,26), (2,45), (4,60), (5,85)]  # what combinations of SSP (first) and RCP (second)?
nensemble = 1000                      # how many ensemble members for the Monte Carlo?
#surgeoptions = [0,1,2]                # which surge data sets to use (0 = original CIAM/DINAS-COAST; 1 = GTSR-corrected D-C; 2 = GTSR nearest data points)
#TESTING:
#ssp_rcp_scenarios = [(5,85)]  # what combinations of SSP (first) and RCP (second)?
#ssp_rcp_scenarios = [(1,26), (2,45), (4,60), (5,85)]  # what combinations of SSP (first) and RCP (second)?
#nensemble = 10                      # how many ensemble members for the Monte Carlo?
surgeoptions = [0]                # which surge data sets to use (0 = original CIAM/DINAS-COAST; 1 = GTSR-corrected D-C; 2 = GTSR nearest data points)

# Now, we actually do the simulations
for surgeoption in surgeoptions
    println("Surge option:",surgeoption,"...")

    for (ssp, rcp) in ssp_rcp_scenarios

        println("Running SSP",ssp, "-RCP",rcp,"...")

        # write the init file
        init_settings = Dict(
            :init_filename   => "MCdriver_init.csv",
            :subset          => false,
            :ssp             => ssp_files[ssp],
            :ssp_simplified  => ssp,
            :rcp             => rcp,
            :b => "lsl_rcp85_p50.csv" # won't matter, overwritten when doing the SLR sampling anyway
        )
        init_file = joinpath(outputdir,init_settings[:init_filename])
        textheader="run_name,lslr,subset,ssp,ssp_simplified\n"
        textstr = "base,$(init_settings[:b]),$(init_settings[:subset]),$(init_settings[:ssp]),$(init_settings[:ssp_simplified])"
        txtfile=open(init_file,"w") do io
            write(io,textheader)
            write(io,textstr)
        end

        # Define input trial parameters: BRICK model, number of trials (n), min and max BRICK percentile,
        # start and end year, and timestep. Resetting each time through the loop for consistency.
        trial_params = Dict(
            :brickfile  => brickfile,
            :n      => nensemble,
            :high   => 100,
            :low    => 0,
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
            :popval         => popinput,
            :GAMSmatch      => false,
            :surgeoption    => 0,
            :subset         => false
        )

        # vary SLR and CIAM parameters
        println("... vary SLR and vary CIAM ...")
        trial_params[:low] = 0
        trial_params[:high] = 100
        runname = string("SSP",init_settings[:ssp_simplified],"_BRICK",init_settings[:rcp],"_surge",surgeoption,"_varySLR_varyCIAM")
        runTrials(init_settings[:rcp], trial_params, adaptRegime1, outputdir, init_file, vary_slr=true, vary_ciam=true, runname=runname)

        # vary only CIAM parameters
        println("... vary CIAM ...")
        trial_params[:low] = 50
        trial_params[:high] = 50
        runname = string("SSP",init_settings[:ssp_simplified],"_BRICK",init_settings[:rcp],"_surge",surgeoption,"_varyCIAM")
        runTrials(init_settings[:rcp], trial_params, adaptRegime1, outputdir, init_file, vary_slr=false, vary_ciam=true, runname=runname)

        # vary only BRICK parameters
        println("... vary SLR ...")
        trial_params[:low] = 0
        trial_params[:high] = 100
        runname = string("SSP",init_settings[:ssp_simplified],"_BRICK",init_settings[:rcp],"_surge",surgeoption,"_varySLR")
        runTrials(init_settings[:rcp], trial_params, adaptRegime1, outputdir, init_file, vary_slr=true, vary_ciam=false, runname=runname)

        # only 5th percentile of SLR, with CIAM defaults
        println("... 5th percentile ...")
        prctile = 5
        trial_params[:low] = prctile
        trial_params[:high] = prctile
        trial_params[:n] = 1
        runname = string("SSP",init_settings[:ssp_simplified],"_BRICK",init_settings[:rcp],"_surge",surgeoption,"_p",prctile)
        runTrials(init_settings[:rcp], trial_params, adaptRegime1, outputdir, init_file, vary_slr=false, vary_ciam=false, runname=runname)

        # only 50th percentile of SLR, with CIAM defaults
        println("... 50th percentile ...")
        prctile = 50
        trial_params[:low] = prctile
        trial_params[:high] = prctile
        trial_params[:n] = 1
        runname = string("SSP",init_settings[:ssp_simplified],"_BRICK",init_settings[:rcp],"_surge",surgeoption,"_p",prctile)
        runTrials(init_settings[:rcp], trial_params, adaptRegime1, outputdir, init_file, vary_slr=false, vary_ciam=false, runname=runname)

        # only 95th percentile of SLR, with CIAM defaults
        println("... 95th percentile ...")
        prctile = 95
        trial_params[:low] = prctile
        trial_params[:high] = prctile
        trial_params[:n] = 1
        runname = string("SSP",init_settings[:ssp_simplified],"_BRICK",init_settings[:rcp],"_surge",surgeoption,"_p",prctile)
        runTrials(init_settings[:rcp], trial_params, adaptRegime1, outputdir, init_file, vary_slr=false, vary_ciam=false, runname=runname)

    end
end

##==============================================================================
## End
##==============================================================================
