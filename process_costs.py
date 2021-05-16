def process_costs_df(dfC, dfO, dfN, tmax=False):
    # set up
    ntime = len(dfO.time.unique())
    levs = {"RetreatCost" : [1,10,100,1000,10000], "ProtectCost" : [10,100,1000,10000], "NoAdaptCost" : [0]}
    
    # figure out which segments pursue which adaptation options, and at what levels of protection
    actions = {}
    for choice in levs.keys():
        for lev in levs[choice]:
            if lev > 0:
                actions[choice+str(lev)] = list(dfO.loc[(dfO.time==1)&(dfO.variable==choice)&(dfO.level==lev),"segments"])
            else:
                actions[choice+str(lev)] = list(dfO.loc[(dfO.time==1)&(dfO.variable==choice),"segments"])
    retreat_segs, protect_segs = [], []
    for lev in levs["RetreatCost"]:
        retreat_segs += actions["RetreatCost"+str(lev)]
    for lev in levs["ProtectCost"]:
        protect_segs += actions["ProtectCost"+str(lev)]
        
    # tally up costs associated with each action
    retreat_costs_N = [0]*ntime
    protect_costs_N = [0]*ntime
    inundation_costs_N = [0]*ntime
    wetland_costs_N = [0]*ntime
    flood_costs_N = [0]*ntime
    
    for t in range(1,ntime+1):
        dfNsub = dfN.loc[(dfN.time==t)] # subset to speed the loop up
        # retreat
        for lev in levs["RetreatCost"]:
            pname = "RetreatCost"+str(lev)
            retreat_costs_N[t-1] += dfNsub.loc[(dfNsub.segments.isin(actions[pname])) & (dfNsub.variable=="RelocateRetreat") & (dfNsub.level==lev), "value"].sum()
        retreat_costs_N[t-1] += dfNsub.loc[(dfNsub.segments.isin(actions["NoAdaptCost0"])) & (dfNsub.variable=="RelocateNoAdapt"), "value"].sum()
        # protect
        for lev in levs["ProtectCost"]:
            pname = "ProtectCost"+str(lev)
            protect_costs_N[t-1] += dfNsub.loc[(dfNsub.segments.isin(actions[pname])) & (dfNsub.variable=="Construct") & (dfNsub.level==lev), "value"].sum()
        # inundation (Flood)
        for lev in levs["RetreatCost"]:
            pname = "RetreatCost"+str(lev)
            inundation_costs_N[t-1] += dfNsub.loc[(dfNsub.segments.isin(actions[pname])) & (dfNsub.variable=="FloodRetreat") & (dfNsub.level==lev), "value"].sum()
        inundation_costs_N[t-1] += dfNsub.loc[(dfNsub.segments.isin(actions["NoAdaptCost0"])) & (dfNsub.variable=="FloodNoAdapt"), "value"].sum()
        # wetland
        for lev in levs["RetreatCost"]:
            pname = "RetreatCost"+str(lev)
            wetland_costs_N[t-1] += dfNsub.loc[(dfNsub.segments.isin(actions[pname])) & (dfNsub.variable=="WetlandRetreat"), "value"].sum()
        for lev in levs["ProtectCost"]:
            pname = "ProtectCost"+str(lev)
            wetland_costs_N[t-1] += dfNsub.loc[(dfNsub.segments.isin(actions[pname])) & (dfNsub.variable=="WetlandProtect"), "value"].sum()
        wetland_costs_N[t-1] += dfNsub.loc[(dfNsub.segments.isin(actions["NoAdaptCost0"])) & (dfNsub.variable=="WetlandNoAdapt"), "value"].sum()
        # flooding (Storm)
        for lev in levs["RetreatCost"]:
            pname = "RetreatCost"+str(lev)
            flood_costs_N[t-1] += dfNsub.loc[(dfNsub.segments.isin(actions[pname])) & (dfNsub.level==lev) & (dfNsub.variable.isin(["StormCapitalRetreat","StormPopRetreat"])), "value"].sum()
        for lev in levs["ProtectCost"]:
            pname = "ProtectCost"+str(lev)
            flood_costs_N[t-1] += dfNsub.loc[(dfNsub.segments.isin(actions[pname])) & (dfNsub.level==lev) & (dfNsub.variable.isin(["StormCapitalProtect","StormPopProtect"])), "value"].sum()
        flood_costs_N[t-1] += dfNsub.loc[(dfNsub.segments.isin(actions["NoAdaptCost0"])) & (dfNsub.variable.isin(["StormCapitalNoAdapt","StormPopNoAdapt"])), "value"].sum()

    if not tmax:
        tmax = ntime
        
    dfNew = pd.DataFrame()
    dfNew["time"] = list(range(1,ntime+1))
    dfNew["NoAdapt"] = list(dfC.loc[(dfC.variable=="NoAdaptCost"), "value"])
    dfNew["Optimal"] = list(dfC.loc[(dfC.variable=="OptimalCost"), "value"])
    dfNew["FloodNoAdapt"] = list(dfC.loc[(dfC.variable=="FloodNoAdapt"), "value"])
    dfNew["WetlandNoAdapt"] = list(dfC.loc[(dfC.variable=="WetlandNoAdapt"), "value"])
    dfNew["RelocateNoAdapt"] = list(dfC.loc[(dfC.variable=="RelocateNoAdapt"), "value"])
    dfNew["StormCapitalNoAdapt"] = list(dfC.loc[(dfC.variable=="StormCapitalNoAdapt"), "value"])
    dfNew["StormPopNoAdapt"] = list(dfC.loc[(dfC.variable=="StormPopNoAdapt"), "value"])
    dfNew["RetreatOptimal"] = retreat_costs_N
    dfNew["ProtectOptimal"] = protect_costs_N
    dfNew["InundationOptimal"] = inundation_costs_N
    dfNew["WetlandOptimal"] = wetland_costs_N
    dfNew["FloodOptimal"] = flood_costs_N
    dfNew = dfNew.loc[(dfNew.time <= tmax)]
    
    return dfNew