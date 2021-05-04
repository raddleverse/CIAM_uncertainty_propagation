# baseline comparisons against the GAMS results

using Mimi, MimiCIAM, Query, RData, StatsBase, CSV, DataFrames, NetCDF

include("brickLSL.jl")

cd("/Users/aewsma/codes/CIAM_adaptation_regimes/ciam-code/src")

##==============================================================================
## baseline cases
noRetreat = false
allowMaintain = false
fixed = true
subset = false
#subset = "segments_test.csv"
#subset = "sub1names.csv"

## no SLR case and RCP8.5 SLR cases (median and 5th/95th percentiles)
#for b in ["lsl_rcp85_p50.csv"]#"lsl_rcp0_p50.csv","lsl_rcp85_p50.csv","lsl_rcp85_p5.csv","lsl_rcp85_p95.csv"]
pop = 0
SSP = 0
SSP_simp = 2 # won't matter, just to get defaults for the popdens_seg_jones and _merkens arrays.
runname="lslOld_pop0"
b= "lsl_rcp85_p50.csv"
textheader="run_name,lslr,subset,ssp,ssp_simplified\n"
textstr = "base,$(b),$(subset),$(SSP[1]),$(SSP_simp)"
txtfile=open("../data/batch/init.txt","w") do io
    write(io,textheader)
    write(io,textstr)
end
m = MimiCIAM.get_model(t=20,initfile="../data/batch/init.txt",fixed=fixed,noRetreat=noRetreat,allowMaintain=allowMaintain)
run(m)
MimiCIAM.write_ciam(m; runname=runname, sumsegs="global", varnames=false, tag="test")
MimiCIAM.write_ciam(m; runname=runname, sumsegs="seg", varnames=false, tag="test")
MimiCIAM.write_optimal_costs(m; runname=runname)
    # need the breakdown of optimal costs since they're different for each segment
    #  (different adaptation strategies)
#    if b in ["lsl_rcp0_p50.csv","lsl_rcp85_p50.csv"]
        #MimiCIAM.write_optimal_costs(m; runname=runname)
#    end
#end
##==============================================================================


##==============================================================================
## changing things

## baseline+updated SLR
rcp = 85
segIDs = false
brickfile = "../data/lslr/BRICK_projections.RData"

# end_year here only for picking median SLR ensemble member
lsl = brick_lsl(rcp,segIDs,brickfile,1,50,50,2010,2150,10,false)
lslr=lsl[1]
gmsl=lsl[2]
ensInds=lsl[3] # Indices of original BRICK array
update_param!(m,:time,15)
update_param!(m,:lslr,lslr[1,:,:])
run(m)
MimiCIAM.write_ciam(m; runname=runname, sumsegs="global", varnames=false, tag="test")
MimiCIAM.write_ciam(m; runname=runname, sumsegs="seg", varnames=false, tag="test")
MimiCIAM.write_optimal_costs(m; runname=runname)


## baseline+updated GDP/POP via SSP5 (and 2, for SOM)

## baseline+updated SLR+updated GPD/POP via SSP5 (and 2, for SOM)


##==============================================================================


##==============================================================================
## End
##==============================================================================
