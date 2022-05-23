##==============================================================================
## baseline_comparisons.jl
## baseline comparisons against the GAMS results
##
## Tony Wong (aewsma@rit.edu)
##==============================================================================

using Mimi
using MimiCIAM
using CSV
using RData
using DataFrames

include("brickLSL.jl")

# Note that this file uses some built in MimiCIAM functions which you may want
# to copy with a new function name and modify if they do not suit your needs. They are
# displayed below for convenience, but also found in the MimiCIAM package utils.jl

# function write_init_file(run_name::String, outputdir::String, init_settings::Dict)
#     textheader="run_name,lslr,subset,ssp,ssp_simplified\n"
#     textstr = "$(run_name),$(init_settings[:lslrfile]),$(init_settings[:subset]),$(init_settings[:ssp]),$(init_settings[:ssp_simplified])"
#     txtfile = open(joinpath(outputdir, init_settings[:init_filename]),"w") do io
#         write(io,textheader)
#         write(io,textstr)
#     end
# end

# function write_output_files(m, outputdir::String, run_name::String)
#     println("Writing out ciam `subsegs = seg` file for run $(run_name) to directory $(outputdir)")
#     MimiCIAM.write_ciam(m; outputdir = outputdir, runname = run_name, sumsegs="seg", varnames=false)
#     println("Writing out ciam `subsegs = global` file for run $(run_name) to directory $(outputdir)")
#     MimiCIAM.write_ciam(m; outputdir = outputdir, runname = run_name, sumsegs="global", varnames=false)
#     println("Writing out optimal costs file for run $(run_name) to directory $(outputdir)")
#     MimiCIAM.write_optimal_costs(m; outputdir = outputdir, runname = run_name)
# end

##==============================================================================
## Setup

outputdir = joinpath(@__DIR__, "..", "output", "baseline_comparisons")
isdir(outputdir) || mkpath(outputdir)
#brickfile = joinpath(@__DIR__, "..", "data", "lslr", "BRICK_projections.RData")
brickfile = "https://zenodo.org/record/6461560/files/sneasybrick_projections_csv.zip"

##==============================================================================
## Helper Functions

"""
    get_segIDs(subset)

Return the segment IDs of a subset of segments.
"""
function get_segIDs(subset)
    if subset == false
        return false
    else
        subs = MimiCIAM.load_subset(subset)
        sort!(subs)
        return MimiCIAM.segStr_to_segID(subs)
    end
end

##==============================================================================
## Run Comparisons

##==============================================================================
## ctrl+noConstrFix: This case is run with a modified slrcost component held in
## slrcost_GAMSmatch.jl, which is taken care of in the `get_model` step with the
## GAMS match arg and removes the block that disallows height reductions

run_name = "ctrl+noConstrFix"

init_settings = Dict(
    :init_filename   => string("$run_name", "_init.csv"),
    :lslrfile        => "lsl_rcp85_p50.csv",
    :subset          => false,
    :ssp             => 0,
    :ssp_simplified  => 2 # won't matter, just to get defaults for the popdens_seg_jones and _merkens arrays.
)

model_settings = Dict(
    :fixed          => true,
    :t              => 20,
    :noRetreat      => false,
    :allowMaintain  => false,
    :popinput       => 0,
    :GAMSmatch      => true,
    :surgeoption    => 0
)

MimiCIAM.write_init_file(run_name, outputdir, init_settings)

m = MimiCIAM.get_model(
    initfile        = joinpath(outputdir, init_settings[:init_filename]),
    fixed           = model_settings[:fixed],
    t               = model_settings[:t],
    noRetreat       = model_settings[:noRetreat],
    allowMaintain   = model_settings[:allowMaintain],
    popinput        = model_settings[:popinput],
    GAMSmatch       = model_settings[:GAMSmatch],
    surgeoption     = model_settings[:surgeoption]
)
run(m)

MimiCIAM.write_output_files(m, outputdir, run_name)

##==============================================================================
##  Ctrl Case

run_name = "ctrl"

init_settings = Dict(
    :init_filename   => string("$run_name", "_init.csv"),
    :lslrfile        => "lsl_rcp85_p50.csv",
    :subset          => false,
    :ssp             => 0,
    :ssp_simplified  => 2 # won't matter, just to get defaults for the popdens_seg_jones and _merkens arrays.
)

model_settings = Dict(
    :fixed          => true,
    :t              => 15,
    :noRetreat      => false,
    :allowMaintain  => false,
    :popinput       => 0,
    :GAMSmatch      => false,
    :surgeoption    => 0
)

MimiCIAM.write_init_file(run_name, outputdir, init_settings)

m = MimiCIAM.get_model(
    initfile        = joinpath(outputdir, init_settings[:init_filename]),
    fixed           = model_settings[:fixed],
    t               = model_settings[:t],
    noRetreat       = model_settings[:noRetreat],
    allowMaintain   = model_settings[:allowMaintain],
    popinput        = model_settings[:popinput],
    GAMSmatch       = model_settings[:GAMSmatch],
    surgeoption     = model_settings[:surgeoption]
)
run(m)

MimiCIAM.write_output_files(m, outputdir, run_name)

##==============================================================================
##  baseline+updated GDP/POP via SSP5. but can be any of 1-5

run_name = "ctrl+SSP5"

init_settings = Dict(
    :init_filename   => string("$run_name", "_init.csv"),
    :lslrfile        => "lsl_rcp85_p50.csv",
    :subset          => false,
    :ssp             => "IIASAGDP_SSP5_v9_130219",
    :ssp_simplified  => 5
)

model_settings = Dict(
    :fixed          => true,
    :t              => 15,
    :noRetreat      => false,
    :allowMaintain  => false,
    :popinput       => 0,
    :GAMSmatch      => false,
    :surgeoption    => 0
)

MimiCIAM.write_init_file(run_name, outputdir, init_settings)

m = MimiCIAM.get_model(
    initfile        = joinpath(outputdir, init_settings[:init_filename]),
    fixed           = model_settings[:fixed],
    t               = model_settings[:t],
    noRetreat       = model_settings[:noRetreat],
    allowMaintain   = model_settings[:allowMaintain],
    popinput        = model_settings[:popinput],
    GAMSmatch       = model_settings[:GAMSmatch],
    surgeoption     = model_settings[:surgeoption]
)
run(m)

MimiCIAM.write_output_files(m, outputdir, run_name)



##==============================================================================
##  baseline+updated SLR

for rcp in [26,45,60,85]

    run_name = "ctrl+BRICKLSL$(rcp)"

    init_settings = Dict(
        :init_filename   => string("$run_name", "_init.csv"),
        :lslrfile        => "lsl_rcp85_p50.csv",
        :subset          => false,
        :ssp             => 0,
        :ssp_simplified  => "2" # won't matter, just to get defaults for the popdens_seg_jones and _merkens arrays.
    )

    model_settings = Dict(
        :fixed          => true,
        :t              => 15,
        :noRetreat      => false,
        :allowMaintain  => false,
        :popinput       => 0, # 0 = default old stuff, 1 = Jones and O'Neill (2016), 2 = Merkens et al (2016)
        :GAMSmatch      => false,
        :surgeoption    => 0
    )

    MimiCIAM.write_init_file(run_name, outputdir, init_settings)

    m = MimiCIAM.get_model(
        initfile        = joinpath(outputdir, init_settings[:init_filename]),
        fixed           = model_settings[:fixed],
        t               = model_settings[:t],
        noRetreat       = model_settings[:noRetreat],
        allowMaintain   = model_settings[:allowMaintain],
        popinput        = model_settings[:popinput],
        GAMSmatch       = model_settings[:GAMSmatch],
        surgeoption     = model_settings[:surgeoption]
    )

    segIDs = get_segIDs(init_settings[:subset])
    lsl = brick_lsl(rcp, segIDs, brickfile, 1, 50, 50, 2010, 2150, 10, false) # end_year here only for picking median SLR ensemble member
    lslr = lsl[1]
    gmsl = lsl[2] # unused
    ensInds = lsl[3] # Indices of original BRICK array (unused)
    update_param!(m, :slrcost, :lslr, lslr[1,:,:])

    run(m)

    MimiCIAM.write_output_files(m, outputdir, run_name)

end

##==============================================================================



##==============================================================================
## Loop over the 4 baseline scenarios, percentiles, and surge options

ssp_files = Dict(1 => "IIASAGDP_SSP1_v9_130219",
                 2 => "IIASAGDP_SSP2_v9_130219",
                 3 => "IIASAGDP_SSP3_v9_130219",
                 4 => "IIASAGDP_SSP4_v9_130219",
                 5 => "IIASAGDP_SSP5_v9_130219")
ssp_rcp_scenarios = [(1,26), (2,45), (4,60), (5,85)]  # what combinations of SSP (first) and RCP (second)?
surgeoptions = [0,1,2]                # which surge data sets to use (0 = original CIAM/DINAS-COAST; 1 = GTSR-corrected D-C; 2 = GTSR nearest data points)
prctiles = [50,5,95]

for surgeoption in surgeoptions
    for (ssp,rcp) in ssp_rcp_scenarios
        for prctile in prctiles

            if prctile == 50
                run_name = "ctrl+SSP$(ssp)+BRICKLSL$(rcp)+surge$(surgeoption)"
            else
                run_name = "ctrl+SSP$(ssp)+BRICKLSL$(rcp)+surge$(surgeoption)_p$(prctile)"
            end

            init_settings = Dict(
                :init_filename   => string("$run_name", "_init.csv"),
                :lslrfile        => "lsl_rcp85_p50.csv", # placeholder - will be overwritten with BRICK LSLR below with `update_param`
                :subset          => false,
                :ssp             => ssp_files[ssp],
                :ssp_simplified  => ssp
            )

            model_settings = Dict(
                :fixed          => true,
                :t              => 15,
                :noRetreat      => false,
                :allowMaintain  => false,
                :popinput       => 0, # 0 = default old stuff, 1 = Jones and O'Neill (2016), 2 = Merkens et al (2016)
                :GAMSmatch      => false,
                :surgeoption    => surgeoption
            )

            MimiCIAM.write_init_file(run_name, outputdir, init_settings)

            m = MimiCIAM.get_model(
                initfile        = joinpath(outputdir, init_settings[:init_filename]),
                fixed           = model_settings[:fixed],
                t               = model_settings[:t],
                noRetreat       = model_settings[:noRetreat],
                allowMaintain   = model_settings[:allowMaintain],
                popinput        = model_settings[:popinput],
                GAMSmatch       = model_settings[:GAMSmatch],
                surgeoption     = model_settings[:surgeoption]
            )

            segIDs = get_segIDs(init_settings[:subset])
            # rcp set in loop now
            lsl = brick_lsl(rcp,segIDs,brickfile,1,prctile,prctile,2010,2150,10,false) # end_year here only for picking median SLR ensemble member
            lslr=lsl[1]
            gmsl=lsl[2]
            ensInds=lsl[3] # Indices of original BRICK array
            update_param!(m, :slrcost, :lslr, lslr[1,:,:])

            run(m)

            MimiCIAM.write_output_files(m, outputdir, run_name)

        end
    end
end

##==============================================================================
## End
##==============================================================================
