
function runTrials(rcp,trial_params,adaptRegime)
    #outfilepath = joinpath("/Volumes/MASTERS/ciammcs","CIAM$(Dates.format(now(), "yyyy-mm-dd HH-MM-SS"))MC$(trial_params["n"])Reg$(adaptRegime["regimeNum"])")
    outfilepath = "../../output"
    mkpath("$outfilepath/monteCarlo_results")

    # Output Files: Trials, NPV, Global Time Series, Regional Spotlight

    # Load CIAM parameters from file
    if trial_params["subset"][1]==false
        segIDs=false
    else
        subs = MimiCIAM.load_subset(trial_params["subset"])
        sort!(subs)
        segIDs = MimiCIAM.segStr_to_segID(subs)
    end

    # Load BRICK data
    lsl = brick_lsl(rcp,segIDs,trial_params["brickfile"],trial_params["n"],trial_params["low"],trial_params["high"],trial_params["ystart"],trial_params["yend"],trial_params["tstep"],false)
    lslr=lsl[1]
    gmsl=lsl[2]
    ensInds=lsl[3] # Indices of original BRICK array
    num_ens=trial_params["n"]

    m = MimiCIAM.get_model(t=trial_params["t"], initfile="../data/batch/init.txt",
                            fixed=adaptRegime["fixed"],noRetreat=adaptRegime["noRetreat"],
                            allowMaintain=adaptRegime["allowMaintain"], popinput=adaptRegime["popval"])
#    update_param!(m,:popinput,adaptRegime["popval"])

    global outtrials=DataFrame()
    global outts=DataFrame()
    globalNPV = zeros(num_ens)
    rgns=["USA"] # USA Seg IDs of interest
    for i=1:num_ens
        update_param!(m,:lslr,lslr[i,:,:])
        res = run_ciam_mcs(m, 1, 10, outfilepath, false)
        res1 = DataFrame([res.current_data])
        global outtrials = [outtrials;res1]
        ts = MimiCIAM.getTimeSeries(m,i,rgns=false,sumsegs="global")
        global outts = [outts;ts]
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
