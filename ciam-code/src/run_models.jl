## Run CIAM-Kopp and CIAM-BRICK 
#using MimiCIAM
include("MimiCIAM.jl")
BRICKslr = ["BRICKsneasy-lsl_rcp85_p5.csv","BRICKsneasy-lsl_rcp85_p10.csv",
"BRICKsneasy-lsl_rcp85_p50.csv","BRICKsneasy-lsl_rcp85_p90.csv","BRICKsneasy-lsl_rcp85_p95.csv"]
KoppSlr = ["lsl_rcp85_p5.csv","lsl_rcp85_p50.csv","lsl_rcp85_p95.csv"]
SSP=["IIASAGDP_SSP5_v9_130219"]
SSP_simp = "5" # TWmod

### Loop: Each run, modify init.txt with BRICK-SSP pairing, initialize model, write optimal results.
runname="BRICKsneasy-IIASAGDP"  
#textheader="run_name,lslr,subset,ssp\n"
textheader="run_name,lslr,subset,ssp,ssp_simplified\n" # TWmod
for b in BRICKslr
#    textstr = "base,$(b),false,$(SSP[1])"
    textstr = "base,$(b),false,$(SSP[1]),$(SSP_simp)" # TWmod
    txtfile=open("../data/batch/init.txt","w") do io
        write(io,textheader)
        write(io,textstr)
    end
    m = MimiCIAM.get_model()
    run(m)
    MimiCIAM.write_optimal_costs(m,runname=runname)

end
    