##==============================================================================
## run_ciam_mcs.jl
##
## Original code: Catherine Ledna (4 Feb 2021)
## Modified code: Tony Wong (16 May 2021)
##==============================================================================

using Dates
using CSV

"""
Run a Monte Carlo simulation with the CIAM model over its distributional parameters.
- trials: number of trials to run.
- ntsteps: number of timesteps to run
- output_dir: an output directory
- save_trials: whether to generate and save MC trial values to a file
"""

# The file structure created by this process will look as follows:
#
# - top directory output_dir containing a RawResults folder with results directly
#   from the monte carlo runs
#

function run_ciam_mcs(model, output_dir; trials=10000, ntsteps=10, save_trials=true, vary_ciam=true)

    # results output directory
    results_output_dir = joinpath(output_dir, "RawResults")
    isdir(results_output_dir) || mkpath(results_output_dir)

    # trials output file
    trials_output_filename = save_trials ? joinpath(outputdir, "trials.csv") : nothing

    # Get an instance of CIAM's mcs
    mcs = getmcs(vary_ciam)

    # run monte carlo trials
    res = run(mcs, model, trials;
        trials_output_filename = trials_output_filename,
        ntimesteps = ntsteps, results_output_dir = results_output_dir)

    return res

end
