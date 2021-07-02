## Run CIAM-Kopp and CIAM-BRICK 
using MimiCIAM

# pull in helper functions
include("utils.jl")

# create output directory
outputdir = joinpath(@__DIR__, "..", "output", "SLR_runs")
isdir(outputdir) || mkdir(outputdir)

## BRICK RUNS

BRICKslr = [
    "BRICKsneasy-lsl_rcp85_p5.csv",
    "BRICKsneasy-lsl_rcp85_p10.csv",
    "BRICKsneasy-lsl_rcp85_p50.csv",
    "BRICKsneasy-lsl_rcp85_p90.csv",
    "BRICKsneasy-lsl_rcp85_p95.csv"
]

# For each run, create a new init.txt with BRICK-SSP pairing, initialize model, write 
# optimal results.
for b in BRICKslr

    run_name = "BRICKsneasy-IIASAGDP_$(split(b, ".")[1])"

    init_settings = Dict(
        :init_filename   => string("$run_name", "_init.csv"),
        :lslrfile        => b,
        :subset          => false,
        :ssp             => "IIASAGDP_SSP5_v9_130219",
        :ssp_simplified  => 5
    )

    MimiCIAM.write_init_file(run_name, outputdir, init_settings)
    m = MimiCIAM.get_model(initfile = joinpath(outputdir, init_settings[:init_filename]))
    run(m)

    MimiCIAM.write_optimal_costs(m; outputdir = outputdir, runname = run_name)

end
    
## KOPP RUNS

KoppSlr = [
    "lsl_rcp0_p50.csv",
    "lsl_rcp45_p50.csv",
    "lsl_rcp85_p5.csv",
    "lsl_rcp85_p50.csv",
    "lsl_rcp85_p95.csv"
]

# For each run, create a new init.txt with BRICK-SSP pairing, initialize model, write 
# optimal results.
for k in KoppSlr

    run_name = "BRICKsneasy-IIASAGDP_$(split(k, ".")[1])"

    init_settings = Dict(
        :init_filename   => string("$run_name", "_init.csv"),
        :lslrfile        => k,
        :subset          => false,
        :ssp             => "IIASAGDP_SSP5_v9_130219",
        :ssp_simplified  => 5
    )

    MimiCIAM.write_init_file(run_name, outputdir, init_settings)
    m = MimiCIAM.get_model(initfile = joinpath(outputdir, init_settings[:init_filename]))
    run(m)

    MimiCIAM.write_optimal_costs(m; outputdir = outputdir, runname = run_name)

end
