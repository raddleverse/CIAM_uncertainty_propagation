
using Mimi
include("src/montecarlo/defmcs.jl")
include("src/montecarlo/run_ciam_mcs.jl")
include("src/brickLSL.jl")
include("src/ciamhelper.jl")
include("src/processResults.jl")

# Define input trial parameters: BRICK model, number of trials (n), min and max BRICK percentile,
#   start and end year, timestep and number of timesteps 
trial_params=Dict(
    "brickfile" => "/Users/catherineledna/Desktop/BRICK/output_model/BRICK-fastdyn_physical_uniform_26May2020.nc",
    "n" => 2880,
    "high" => 100,
    "low"=>0,
    "ystart"=>2010,
    "yend"=>2100,
    "tstep"=>10,
    "t"=>10
)

adaptRegime1=Dict(
    "allowMaintain"=>false,
    "noRetreat"=>false,
    "popval"=>1,
    "fixed"=>false,
    "regimeNum"=>1
)

adaptRegime2=Dict(
    "allowMaintain"=>false,
    "noRetreat"=>true,
    "popval"=>1,
    "fixed"=>false,
    "regimeNum"=>2
)

adaptRegime3=Dict(
    "allowMaintain"=>false,
    "noRetreat"=>false,
    "popval"=>1,
    "fixed"=>true,
    "regimeNum"=>3
)

adaptRegime4=Dict(
    "allowMaintain"=>false,
    "noRetreat"=>true,
    "popval"=>1,
    "fixed"=>true,
    "regimeNum"=>4
)

function runTrials(rcp,trial_params,adaptRegime)
    outfilepath = joinpath("/Volumes/MASTERS/ciammcs","CIAM$(Dates.format(now(), "yyyy-mm-dd HH-MM-SS"))MC$(trial_params["n"])Reg$(adaptRegime["regimeNum"])")
    mkpath("$outfilepath/results")

    # Output Files: Trials, NPV, Global Time Series, Regional Spotlight  

    # Load CIAM parameters from file 
    ciamparams = MimiCIAM.init()
    if ciamparams["subset"][1]==false
        segIDs=false
    else
        subs = MimiCIAM.load_subset(ciamparams["subset"])
        sort!(subs)
        segIDs = MimiCIAM.segStr_to_segID(subs)
    end

    # Load BRICK data 
    lsl = brick_lsl(rcp,segIDs,trial_params["brickfile"],trial_params["n"],trial_params["low"],trial_params["high"],trial_params["ystart"],trial_params["yend"],trial_params["tstep"],false)
    lslr=lsl[1]
    gmsl=lsl[2]
    ensInds=lsl[3] # Indices of original BRICK array
    num_ens=trial_params["n"]
    
    m = MimiCIAM.get_model(t=trial_params["t"],fixed=adaptRegime["fixed"],noRetreat=adaptRegime["noRetreat"],allowMaintain=adaptRegime["allowMaintain"])
    update_param!(m,:popinput,adaptRegime["popval"])

    global outtrials=DataFrame() 
    global outts=DataFrame()
    globalNPV = zeros(num_ens)
    rgns=["USA"] # USA Seg IDs of interest 
    for i=1:num_ens
        update_param!(m,:lslr,lslr[i,:,:])
        res=run_ciam_mcs(m,1,10,outfilepath, false)

        res1=DataFrame([res.current_data])
        global outtrials=[outtrials;res1]
        ts= getTimeSeries(m,i,rgns=rgns)
        global outts=[outts;ts]
        globalNPV[i] = m[:slrcost,:NPVOptimalTotal]
    end

    # Write Trials, Global NPV and Time Series 
    outtrialsname="$(outfilepath)/results/trials.csv"
    CSV.write(outtrialsname,outtrials)

    outnpvname="$(outfilepath)/results/globalnpv.csv"
    procGlobalOutput(globalNPV,gmsl,ensInds,trial_params["brickfile"],rcp,adaptRegime["noRetreat"],outnpvname)

    outtsname="$(outfilepath)/results/globalts_$(rcp).csv"
    CSV.write(outtsname,outts)

end

function runTrialsParallel(rcp,lsl,chunk,trial_params,adaptRegime)
    outfilepath = joinpath("/Volumes/MASTERS/ciammcs","CIAM$(Dates.format(now(), "yyyy-mm-dd HH-MM-SS"))MC$(trial_params["n"])Reg$(adaptRegime["regimeNum"])Chunk$(chunk)")
    mkpath("$outfilepath/results")

    # Output Files: Trials, NPV, Global Time Series, Regional Spotlight  

    # Load CIAM parameters from file 
    ciamparams = MimiCIAM.init()
    if ciamparams["subset"][1]==false
        segIDs=false
    else
        subs = MimiCIAM.load_subset(ciamparams["subset"])
        sort!(subs)
        segIDs = MimiCIAM.segStr_to_segID(subs)
    end

    # Load BRICK data 
    #lsl = brick_lsl(rcp,segIDs,trial_params["brickfile"],trial_params["n"],trial_params["low"],trial_params["high"],trial_params["ystart"],trial_params["yend"],trial_params["tstep"],false)
    lslr=lsl[1]
    gmsl=lsl[2]
    ensInds=lsl[3] # Indices of original BRICK array
    num_ens=size(lslr)[1]   #trial_params["n"]
    
    m = MimiCIAM.get_model(t=trial_params["t"],fixed=adaptRegime["fixed"],noRetreat=adaptRegime["noRetreat"],allowMaintain=adaptRegime["allowMaintain"])
    update_param!(m,:popinput,adaptRegime["popval"])

    global outtrials=DataFrame() 
    global outts=DataFrame()
    globalNPV = zeros(num_ens)
    rgns=["USA"] # USA Seg IDs of interest 
    for i=1:num_ens
        update_param!(m,:lslr,lslr[i,:,:])
        res=run_ciam_mcs(m,1,10,outfilepath, false)

        res1=DataFrame([res.current_data])
        global outtrials=[outtrials;res1]
        ts= getTimeSeries(m,i,rgns=rgns)
        global outts=[outts;ts]
        globalNPV[i] = m[:slrcost,:NPVOptimalTotal]

    end

    # Write Trials, Global NPV and Time Series 
    outtrialsname="$(outfilepath)/results/trials.csv"
    CSV.write(outtrialsname,outtrials)

    outnpvname="$(outfilepath)/results/globalnpv.csv"
    procGlobalOutput(globalNPV,gmsl,ensInds,trial_params["brickfile"],rcp,adaptRegime["noRetreat"],outnpvname)

    outtsname="$(outfilepath)/results/globalts_$(rcp).csv"
    CSV.write(outtsname,outts)

end

runTrials(85,trial_params,adaptRegime1)
runTrials(45,trial_params,adaptRegime1)
runTrials(85,trial_params,adaptRegime2)
runTrials(85,trial_params,adaptRegime3)

runTrials(85,trial_params,adaptRegime4)

# Attempt at running parallel processes 
# lsl=brick_lsl(45,segIDs,trial_params["brickfile"],trial_params["n"],trial_params["low"],trial_params["high"],trial_params["ystart"],trial_params["yend"],trial_params["tstep"],false)
# lslr1=lsl[1][1:1440,:,:]
# gmsl1=lsl[2][:,1:1440]
# ensInds1=lsl[3][1:1440]
# lsl1=(lslr1,gmsl1,ensInds1)

# lslr2=lsl[1][1441:2880,:,:]
# gmsl2=lsl[2][:,1441:2880]
# ensInds2=lsl[3][1441:2880]
# lsl2=(lslr2,gmsl2,ensInds2)

# runTrialsParallel(45,lsl1,1,trial_params,adaptRegime1)
# runTrialsParallel(45,lsl2,2,trial_params,adaptRegime1)




