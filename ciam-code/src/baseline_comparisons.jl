##==============================================================================
## baseline_comparisons.jl
## baseline comparisons against the GAMS results
##
## Tony Wong (aewsma@rit.edu)
##==============================================================================

using MimiCIAM
using CSV
using RData

include("brickLSL.jl")

##==============================================================================
## Setup

outputdir = joinpath(@__DIR__, "..", "output", "baseline_comparisons.jl")
isdir(outputdir) || mkdir(outputdir)
brickfile = "/Users/lisarennels/JuliaProjects/CIAMPaper/local-data/BRICK_projections.RData"

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
    :GAMSmatch      => true
)

# write files
MimiCIAM.write_init_file(run_name, outputdir, init_settings)

m = MimiCIAM.get_model(
    initfile        = joinpath(outputdir, init_settings[:init_filename]),
    fixed           = model_settings[:fixed],
    t               = model_settings[:t], 
    noRetreat       = model_settings[:noRetreat],
    allowMaintain   = model_settings[:allowMaintain],
    popinput        = model_settings[:popinput],
    GAMSmatch       = model_settings[:GAMSmatch]
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
    :GAMSmatch      => false
)

MimiCIAM.write_init_file(run_name, outputdir, init_settings)

m = MimiCIAM.get_model(
    initfile        = joinpath(outputdir, init_settings[:init_filename]),
    fixed           = model_settings[:fixed],
    t               = model_settings[:t], 
    noRetreat       = model_settings[:noRetreat],
    allowMaintain   = model_settings[:allowMaintain],
    popinput        = model_settings[:popinput],
    GAMSmatch       = model_settings[:GAMSmatch]
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
    :GAMSmatch      => false
)

MimiCIAM.write_init_file(run_name, outputdir, init_settings)

m = MimiCIAM.get_model(
    initfile        = joinpath(outputdir, init_settings[:init_filename]),
    fixed           = model_settings[:fixed],
    t               = model_settings[:t], 
    noRetreat       = model_settings[:noRetreat],
    allowMaintain   = model_settings[:allowMaintain],
    popinput        = model_settings[:popinput],
    GAMSmatch       = model_settings[:GAMSmatch]
)
run(m)

MimiCIAM.write_output_files(m, outputdir, run_name)

##==============================================================================
##  baseline+updated population density from Jones and O'Neill (2016)

run_name = "ctrl+SSP5+popJones"

init_settings = Dict(
    :init_filename   => string("$run_name", "_init.csv"),
    :lslrfile        => "lsl_rcp85_p50.csv",
    :subset          => false,
    :ssp             => "IIASAGDP_SSP5_v9_130219",
    :ssp_simplified  => "5" # based on SSPs, so need to use an SSP
)

model_settings = Dict(
    :fixed          => true,
    :t              => 15,
    :noRetreat      => false,
    :allowMaintain  => false,
    :popinput       => 1, # 0 = default old stuff, 1 = Jones and O'Neill (2016), 2 = Merkens et al (2016)
    :GAMSmatch      => false
)

MimiCIAM.write_init_file(run_name, outputdir, init_settings)

m = MimiCIAM.get_model(
    initfile        = joinpath(outputdir, init_settings[:init_filename]),
    fixed           = model_settings[:fixed],
    t               = model_settings[:t], 
    noRetreat       = model_settings[:noRetreat],
    allowMaintain   = model_settings[:allowMaintain],
    popinput        = model_settings[:popinput],
    GAMSmatch       = model_settings[:GAMSmatch]
)
run(m)

MimiCIAM.write_output_files(m, outputdir, run_name)

##==============================================================================
##  baseline+updated SLR RCP8.5

run_name = "ctrl+BRICKLSL85"

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
    :GAMSmatch      => false
)

MimiCIAM.write_init_file(run_name, outputdir, init_settings)

m = MimiCIAM.get_model(
    initfile        = joinpath(outputdir, init_settings[:init_filename]),
    fixed           = model_settings[:fixed],
    t               = model_settings[:t], 
    noRetreat       = model_settings[:noRetreat],
    allowMaintain   = model_settings[:allowMaintain],
    popinput        = model_settings[:popinput],
    GAMSmatch       = model_settings[:GAMSmatch]
)

segIDs = get_segIDs(init_settings[:subset])
rcp = 85
lsl = brick_lsl(rcp,segIDs,brickfile,1,50,50,2010,2150,10,false) # end_year here only for picking median SLR ensemble member
lslr = lsl[1] 
gmsl = lsl[2] # unused
ensInds = lsl[3] # Indices of original BRICK array (unused)
update_param!(m,:lslr,lslr[1,:,:])

run(m)

MimiCIAM.write_output_files(m, outputdir, run_name)

##==============================================================================
##  baseline+updated SLR RCP2.6

run_name = "ctrl+BRICKLSL26"

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
    :GAMSmatch      => false
)

MimiCIAM.write_init_file(run_name, outputdir, init_settings)

m = MimiCIAM.get_model(
    initfile        = joinpath(outputdir, init_settings[:init_filename]),
    fixed           = model_settings[:fixed],
    t               = model_settings[:t], 
    noRetreat       = model_settings[:noRetreat],
    allowMaintain   = model_settings[:allowMaintain],
    popinput        = model_settings[:popinput],
    GAMSmatch       = model_settings[:GAMSmatch]
)

segIDs = get_segIDs(init_settings[:subset])
rcp = 26
lsl = brick_lsl(rcp,segIDs,brickfile,1,50,50,2010,2150,10,false) # end_year here only for picking median SLR ensemble member
lslr=lsl[1]
gmsl=lsl[2] # unused
ensInds=lsl[3] # Indices of original BRICK array (unused)
update_param!(m,:lslr,lslr[1,:,:])

run(m)

MimiCIAM.write_output_files(m, outputdir, run_name)

##==============================================================================
##  baseline+updated SLR (RCP8.5)+updated GPD/POP via SSP5+updated population density from Jones and O'Neill (2016)

run_name = "ctrl+SSP5+popJones+BRICKLSL85"

init_settings = Dict(
    :init_filename   => string("$run_name", "_init.csv"),
    :lslrfile        => "lsl_rcp85_p50.csv",
    :subset          => false,
    :ssp             => "IIASAGDP_SSP5_v9_130219",
    :ssp_simplified  => "5"  #SSP = "IIASAGDP_SSP5_v9_130219"
)

model_settings = Dict(
    :fixed          => true,
    :t              => 15,
    :noRetreat      => false,
    :allowMaintain  => false,
    :popinput       => 1, # 0 = default old stuff, 1 = Jones and O'Neill (2016), 2 = Merkens et al (2016)
    :GAMSmatch      => false
)

MimiCIAM.write_init_file(run_name, outputdir, init_settings)

m = MimiCIAM.get_model(
    initfile        = joinpath(outputdir, init_settings[:init_filename]),
    fixed           = model_settings[:fixed],
    t               = model_settings[:t], 
    noRetreat       = model_settings[:noRetreat],
    allowMaintain   = model_settings[:allowMaintain],
    popinput        = model_settings[:popinput],
    GAMSmatch       = model_settings[:GAMSmatch]
)

segIDs = get_segIDs(init_settings[:subset])
rcp = 85
lsl = brick_lsl(rcp,segIDs,brickfile,1,50,50,2010,2150,10,false) # end_year here only for picking median SLR ensemble member
lslr=lsl[1]
gmsl=lsl[2]
ensInds=lsl[3] # Indices of original BRICK array
update_param!(m,:lslr,lslr[1,:,:])

run(m)

MimiCIAM.write_output_files(m, outputdir, run_name)

##==============================================================================
##  "ctrl+SSP5+popJones+BRICKLSL85_p05" same as above, but ensemble member giving 
##  the 5th percentile of GMSL used

run_name = "ctrl+SSP5+popJones+BRICKLSL85_p05"

init_settings = Dict(
    :init_filename   => string("$run_name", "_init.csv"),
    :lslrfile        => "lsl_rcp85_p50.csv",
    :subset          => false,
    :ssp             => "IIASAGDP_SSP5_v9_130219",
    :ssp_simplified  => "5"  #SSP = "IIASAGDP_SSP5_v9_130219"
)

model_settings = Dict(
    :fixed          => true,
    :t              => 15,
    :noRetreat      => false,
    :allowMaintain  => false,
    :popinput       => 1, # 0 = default old stuff, 1 = Jones and O'Neill (2016), 2 = Merkens et al (2016)
    :GAMSmatch      => false
)

MimiCIAM.write_init_file(run_name, outputdir, init_settings)

m = MimiCIAM.get_model(
    initfile        = joinpath(outputdir, init_settings[:init_filename]),
    fixed           = model_settings[:fixed],
    t               = model_settings[:t], 
    noRetreat       = model_settings[:noRetreat],
    allowMaintain   = model_settings[:allowMaintain],
    popinput        = model_settings[:popinput],
    GAMSmatch       = model_settings[:GAMSmatch]
)

segIDs = get_segIDs(init_settings[:subset])
rcp = 85
lsl = brick_lsl(rcp,segIDs,brickfile,1,5,5,2010,2150,10,false) # end_year here only for picking median SLR ensemble member
lslr=lsl[1]
gmsl=lsl[2]
ensInds=lsl[3] # Indices of original BRICK array
update_param!(m, :lslr, lslr[1,:,:])

run(m)

MimiCIAM.write_output_files(m, outputdir, run_name)

##==============================================================================
##  ctrl+SSP5+popJones+BRICKLSL85_p95

run_name = "ctrl+SSP5+popJones+BRICKLSL85_p95"

init_settings = Dict(
    :init_filename   => string("$run_name", "_init.csv"),
    :lslrfile        => "lsl_rcp85_p50.csv",
    :subset          => false,
    :ssp             => "IIASAGDP_SSP5_v9_130219",
    :ssp_simplified  => "5"  #SSP = "IIASAGDP_SSP5_v9_130219"
)

model_settings = Dict(
    :fixed          => true,
    :t              => 15,
    :noRetreat      => false,
    :allowMaintain  => false,
    :popinput       => 1, # 0 = default old stuff, 1 = Jones and O'Neill (2016), 2 = Merkens et al (2016)
    :GAMSmatch      => false
)

MimiCIAM.write_init_file(run_name, outputdir, init_settings)

m = MimiCIAM.get_model(
    initfile        = joinpath(outputdir, init_settings[:init_filename]),
    fixed           = model_settings[:fixed],
    t               = model_settings[:t], 
    noRetreat       = model_settings[:noRetreat],
    allowMaintain   = model_settings[:allowMaintain],
    popinput        = model_settings[:popinput],
    GAMSmatch       = model_settings[:GAMSmatch]
)

segIDs = get_segIDs(init_settings[:subset])
rcp = 85
lsl = brick_lsl(rcp,segIDs,brickfile,1,95,95,2010,2150,10,false) # end_year here only for picking median SLR ensemble member
lslr=lsl[1]
gmsl=lsl[2]
ensInds=lsl[3] # Indices of original BRICK array
update_param!(m, :lslr, lslr[1,:,:])

run(m)

MimiCIAM.write_output_files(m, outputdir, run_name)

##==============================================================================
## End
##==============================================================================
