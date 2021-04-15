# baseline comparisons against the GAMS results

using MimiCIAM

cd("/Users/aewsma/codes/CIAM_adaptation_regimes/ciam-code/src")

# baseline case
pop = 0
noRetreat = false
allowMaintain = false
fixed = true
SSP = 0
SSP_simp = 2 # won't matter, just to get defaults for the popdens_seg_jones and _merkens arrays.
subset = false
b = "lsl_rcp0_p50.csv"
textheader="run_name,lslr,subset,ssp,ssp_simplified\n"
textstr = "base,$(b),$(subset),$(SSP[1]),$(SSP_simp)"
txtfile=open("../data/batch/init.txt","w") do io
    write(io,textheader)
    write(io,textstr)
end
#m = MimiCIAM.get_model(t=10,fixed=fixed,noRetreat=noRetreat,allowMaintain=allowMaintain)
m = MimiCIAM.get_model(t=10,initfile="../data/batch/init.txt",fixed=fixed,noRetreat=noRetreat,allowMaintain=allowMaintain)
run(m)

runname="lsl0_baseline"
write_optimal_costs(m;runname=runname)
