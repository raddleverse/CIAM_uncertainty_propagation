# baseline comparisons against the GAMS results

using MimiCIAM

cd("/Users/aewsma/codes/CIAM_adaptation_regimes/ciam-code/src")

##==============================================================================
## baseline cases
pop = 0
noRetreat = false
allowMaintain = false
fixed = true
SSP = 0
SSP_simp = 2 # won't matter, just to get defaults for the popdens_seg_jones and _merkens arrays.
#subset = false
subset = "segments_test.csv"
#subset = "sub1names.csv"
runname="lslOld_pop0"

## no SLR case and RCP8.5 SLR cases (median and 5th/95th percentiles)
for b in ["lsl_rcp85_p50.csv"]#"lsl_rcp0_p50.csv","lsl_rcp85_p50.csv","lsl_rcp85_p5.csv","lsl_rcp85_p95.csv"]
    textheader="run_name,lslr,subset,ssp,ssp_simplified\n"
    textstr = "base,$(b),$(subset),$(SSP[1]),$(SSP_simp)"
    txtfile=open("../data/batch/init.txt","w") do io
        write(io,textheader)
        write(io,textstr)
    end
    m = MimiCIAM.get_model(t=20,initfile="../data/batch/init.txt",fixed=fixed,noRetreat=noRetreat,allowMaintain=allowMaintain)
    run(m)
    #MimiCIAM.write_ciam(m; runname=runname, sumsegs="global", varnames=false, tag="test")
    MimiCIAM.write_ciam(m; runname=runname, sumsegs="seg", varnames=false, tag="test")
    # need the breakdown of optimal costs since they're different for each segment
    #  (different adaptation strategies)
#    if b in ["lsl_rcp0_p50.csv","lsl_rcp85_p50.csv"]
        MimiCIAM.write_optimal_costs(m; runname=runname)
#    end
end
##==============================================================================



#b = "BRICKsneasy-lsl_rcp85_p50.csv"

#MimiCIAM.write_optimal_costs(m; runname=runname)
#MimiCIAM.write_optimal_protect_retreat(m; runname=runname)
MimiCIAM.write_ciam(m; runname=runname, sumsegs="global", varnames=false, tag="")
#MimiCIAM.write_ciam(m; runname=runname, sumsegs="seg", varnames=false, tag="")
#MimiCIAM.write_ciam(m; runname=runname, sumsegs="rgn", varnames=false, tag="")
