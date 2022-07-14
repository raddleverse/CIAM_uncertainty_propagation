##==============================================================================
## ciamMonteCarlo.jl
##
## Original code: Catherine Ledna (4 Feb 2021)
## Modified code: Tony Wong (16 May 2021)
##==============================================================================

# The file structure created by this process will look as follows:
#
# - a (unique) top directory created with the code:
#   joinpath(outputdir, runname, "CIAM $(Dates.format(now(), "yyyy-mm-dd HH-MM-SS")) MC$(trial_params[:n])")
#   which will hold a RawResults folder with results directly from the monte carlo
#   runs in run_ciam_mcs.jl and a PostProcessing folder written with the code at
#   the bottom of ciamMonteCarlo.jl

function runTrials(rcp, ssp, trial_params, adaptRegime, outputdir, init_filepath; vary_slr = true, vary_ciam=true, runname="default_run")


    # make an output directory for this SSP-RCP scenario
    scenario = "SSP$(ssp)_BRICK$(rcp)"
    scenario_outputdir = joinpath(outputdir, scenario)
    isdir(scenario_outputdir) || mkpath(scenario_outputdir)

    # within the scenario output directory, make one for this particular set of simulations
    #outputdir = joinpath(scenario_outputdir, runname, "CIAM $(Dates.format(now(), "yyyy-mm-dd HH-MM-SS")) MC$(trial_params[:n])")
    outputdir = joinpath(scenario_outputdir, runname, "CIAM MC$(trial_params[:n])")
    isdir(outputdir) || mkpath(outputdir)

    # Output Files: Trials, NPV, Global Time Series, Regional Spotlight

    # Load CIAM parameters from file
    if adaptRegime[:subset] == false
        segIDs = false
    else
        subs = MimiCIAM.load_subset(adaptRegime[:subset])
        sort!(subs)
        segIDs = MimiCIAM.segStr_to_segID(subs)
    end

    # Load BRICK data
    lsl = brick_lsl(rcp, segIDs, trial_params[:brickfile], trial_params[:n], trial_params[:low],
                        trial_params[:high],trial_params[:ystart], trial_params[:yend],
                        trial_params[:tstep], false)
    if vary_slr
        lslr    =lsl[1]
        gmsl    =lsl[2]
        ensInds =lsl[3] # Indices of original BRICK array

    elseif trial_params[:low] == trial_params[:high]
        lslr = repeat(lsl[1],outer=(trial_params[:n], 1, 1))
        gmsl = repeat(lsl[2], outer = (1,trial_params[:n]))
        ensInds = fill(lsl[3],trial_params[:n]) # only 1 element coming back so `fill` instead of `repeat`

    else
        error("Not varying SLR in Monte Carlo sampling, but the low and high quantiles requested are not equal.")
    end

    num_ens = trial_params[:n]

    m = MimiCIAM.get_model(t = adaptRegime[:t], initfile = init_filepath,
                            fixed = adaptRegime[:fixed], noRetreat = adaptRegime[:noRetreat],
                            allowMaintain = adaptRegime[:allowMaintain], popinput = adaptRegime[:popval],
                            surgeoption = adaptRegime[:surgeoption])
#    update_param!(m,:popinput,adaptRegime[:popval])

    # get the segments and their corresponding World Bank regions
    dfSR = CSV.read("../data/segments_regions_WB.csv", DataFrame)
    dfSR[!,"ids"] = [parse(Int64,replace(i, r"[^0-9]"=> "")) for i in dfSR[!,"segments"]]

    # unique World Bank regions
    wbrgns = unique(dfSR[!,"global region"])

    global outtrials = DataFrame()
    global outts = DataFrame()
    globalNPV = zeros(num_ens)
    regionNPV = zeros(num_ens, length(wbrgns))
    regionNPV1 = zeros(num_ens, length(wbrgns)) # for first time step optimal costs, to subtract off later

    for i = 1:num_ens
        update_param!(m,:slrcost, :lslr,lslr[i,:,:])

        res = run_ciam_mcs(m, outputdir; trials = 1, ntsteps = adaptRegime[:t], save_trials = false, vary_ciam = vary_ciam)
        res1 = DataFrame([res.current_data])
        global outtrials = [outtrials;res1]

        ts = MimiCIAM.getTimeSeries(m,i,rgns=false,sumsegs="global")
        global outts = [outts;ts]
        globalNPV[i] = m[:slrcost,:NPVOptimalTotal]

        # get the segments for each region and aggregate
        for rgn in wbrgns
            segIDs_rgn = filter(:"global region" => ==(rgn), dfSR)[!,"ids"]
            idx_rgn = findall(x -> x in segIDs_rgn, dfSR[!,"ids"])
            col_rgn = findfirst(x->x==rgn, wbrgns)
            regionNPV[i,col_rgn] = sum(m[:slrcost,:NPVOptimal][idx_rgn])
            regionNPV1[i,col_rgn] = sum(m[:slrcost,:OptimalCost][1,idx_rgn]) # regional (World Bank regions) first time step optimal costs
        end
    end

    # get regional NPV as DataFrame for output
    outregionNPV = DataFrame(regionNPV, :auto)
    rename!(outregionNPV,wbrgns)

    # regional (World Bank regions) first time step optimal costs
    outregionNPV1 = DataFrame(regionNPV1, :auto)
    rename!(outregionNPV1,wbrgns)

    # Write Trials, Global NPV and Time Series
    postprocessing_outputdir = joinpath(outputdir, "PostProcessing")
    isdir(postprocessing_outputdir) || mkpath(postprocessing_outputdir)

    outtrialsname = joinpath(postprocessing_outputdir, "trials_$(runname).csv")
    outnpvname= joinpath(postprocessing_outputdir, "globalnpv_$(runname).csv")
    outrgnname = joinpath(postprocessing_outputdir, "regionnpv_$(runname).csv")
    outtsname= joinpath(postprocessing_outputdir, "globalts_$(rcp)_$(runname).csv")

    CSV.write(outtrialsname, outtrials)
    procGlobalOutput(globalNPV,gmsl,ensInds,trial_params[:brickfile],rcp,adaptRegime[:noRetreat],outnpvname)
    CSV.write(outtsname, outts)
    CSV.write(outrgnname, outregionNPV)

    # also write out the regional (World Bank regions) first time step optimal costs
    outrgnname1= joinpath(postprocessing_outputdir, "regionts_$(runname).csv")
    CSV.write(outrgnname1, outregionNPV1)

end

##==============================================================================
## End
##==============================================================================
