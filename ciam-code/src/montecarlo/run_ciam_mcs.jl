##==============================================================================
## run_ciam_mcs.jl
##
## Original code: Catherine Ledna (4 Feb 2021)
## Modified code: Tony Wong (16 May 2021)
##==============================================================================

using Dates
using CSV

"""
Adapted from run_fund_mcs.jl in MimiFUND.
Run a Monte Carlo simulation with the CIAM model over its distributional parameters.
trials: number of trials to run.
ntsteps: number of timesteps to run
output_dir: an output directory
save_trials: whether to generate and save MC trial values to a file
"""

function run_ciam_mcs(model,trials=10000,ntsteps=10,output_dir=nothing, save_trials=true)
    output_dir = output_dir === nothing ? joinpath(@__DIR__, "../../output/ciammcs", "CIAM $(Dates.format(now(), "yyyy-mm-dd HH-MM-SS")) MC$trials") : output_dir
    mkpath("$output_dir/results")

    trials_output_filename = save_trials ?  joinpath("$output_dir/trials.csv") :  nothing

    # Get an instance of CIAM's mcs
    mcs = getmcs()

    # run monte carlo trials
    res = run(mcs, model, trials;
        trials_output_filename = trials_output_filename,
        ntimesteps = ntsteps, results_output_dir = "$output_dir/results")

    return res

end
