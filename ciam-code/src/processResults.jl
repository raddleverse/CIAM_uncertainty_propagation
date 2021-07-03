##==============================================================================
## processResults.jl
##
## Original code: Catherine Ledna (4 Feb 2021)
## Modified code: Tony Wong (16 May 2021)
##==============================================================================

# Functions to process CIAM data and make plots

function procGlobalOutput(glob,gmsl,inds,brickfile,rcp,noRetreat,outfile=false,tstart=2010,tend=2100)

    # Get Global NPV
    npv = glob

    #gmsl_fd = gmsl[10,fd_inds]
    #gmsl_nofd=gmsl[10,nofd_inds]
    gmsl2040 = gmsl[4,:]
    gmsl2060 = gmsl[6,:]
    gmsl2080 = gmsl[8,:]
    gmsl2100 = gmsl[10,:]
    # CAUTION! this is hardwired for 2010-2100

    if noRetreat==false
        lab="With Retreat"
    else
        lab="No Retreat"
    end

    df1 = DataFrame([npv gmsl2100], :auto)
    #df1[:brick]="Fast Dynamics"
    df1[!,:retreat] = fill(lab, size(df1)[1])

    #df2 = DataFrame([npv_nofd gmsl_nofd], :auto)
    #df2[:brick]="No Fast Dynamics"
    #df2[:retreat]=lab

    #outdf = [df1;df2]
    outdf = df1

    rename!(outdf,[:npv,:gmsl2100,:retreat])
    outdf[!,:brickEnsInd] = inds

    if outfile==false
        CSV.write("output/ciammcs/globalNPV_rcp$(rcp)_noRetreat$(noRetreat)_mcs_newPop_maintainFalse.csv",outdf)
    else
        CSV.write(outfile,outdf)
    end

end

function procSegResults(m,seg,lsl,inds,brickfile,rcp,noRetreat,tstart=2010,tend=2100)
    segmap = load_segmap()
    segIDs = m[:slrcost,:segID]
    segNames =segID_to_seg(segIDs,segmap)

    ensInds_fd = findall(x->x in fd_inds,inds)
    ensInds_nofd = findall(x ->!(x in fd_inds),inds)
    indsFd = inds[ensInds_fd]
    indsNoFd = inds[ensInds_nofd]

    npv=seg[1]
    opt2050=seg[2]
    lev2050=seg[3]
    opt2100=seg[4]
    lev2100=seg[5]
  #  lsl2050=lsl45[:,5,:]
  #  lsl2100=lsl45[:,10,:]

    lsl2050fd = lsl2050[indsFd,:]
    lsl2050nofd = lsl2050[indsNoFd,:]
    lsl2100fd = lsl2100[indsFd,:]
    lsl2100nofd = lsl2100[indsNoFd,:]

    gmslFd = gmsl[indsFd]
    gmslNoFd = gmsl[indsNoFd]

    npvFd = npv[indsFd,:]
    npvNoFd = npv[indsNoFd,:]

    opt2050_Fd = opt2050[indsFd,:]
    opt2050_noFd = opt2050[indsNoFd,:]

    lev2050_Fd = lev2050[indsFd,:]
    lev2050_noFd = lev2050[indsNoFd,:]

    opt2100_Fd = opt2100[indsFd,:]
    opt2100_noFd = opt2100[indsNoFd,:]

    lev2100_Fd = lev2100[indsFd,:]
    lev2100_noFd = lev2100[indsNoFd,:]


    lab1=vcat(fill("No Fast Dynamics",length(nofd_inds)))
    lab2 = vcat(fill("Fast Dynamics",length(fd_inds)))

    if noRetreat==false
        lab3=vcat(fill("With Retreat",length(nofd_inds)))
        lab4=vcat(fill("With Retreat",length(fd_inds)))
    else
        lab3=vcat(fill("No Retreat",length(nofd_inds)))
        lab4=vcat(fill("No Retreat",length(fd_inds)))
    end

    dict1=Dict("npvFd"=>npvFd,"opt2050_Fd"=>opt2050_Fd,
        "lev2050_Fd"=>lev2050_Fd,"opt2100_Fd"=>opt2100_Fd,
        "lev2100_Fd"=>lev2100_Fd,"npvNoFd"=>npvNoFd,
        "opt2050_noFd"=>opt2050_noFd,"lev2050_noFd"=>lev2050_noFd,
        "opt2100_noFd"=>opt2100_noFd,"lev2100_noFd"=>lev2100_noFd)
    #dict2=Dict("lsl2050fd"=> lsl2050fd,"lsl2050nofd"=>lsl2050nofd,"lsl2100fd"=>lsl2100fd,
    #"lsl2100nofd"=>lsl2100nofd)

    for k in string.(keys(dict1))
        val = dict1[k]
        if k in ["npvFd","npvNoFd","lsl2050fd","lsl2050nofd","lsl2100fd","lsl2100nofd"]
            mn = mean(val,dims=1)'
            sdev = std(val,dims=1)'
            minval = minimum(val,dims=1)'
            maxval = maximum(val,dims=1)'
            pctl95 = [percentile(val[:,i],95) for i in 1:size(val)[2]]
            pctl5 = [percentile(val[:,i],5) for i in 1:size(val)[2]]
            med = median(val,dims=1)'

            df = DataFrame([segNames mn med sdev minval maxval pctl95 pctl5], :auto)
            names!(df,[:segment,:mean,:median,:stdev,:min,:max,:pct95,:pct5])
            CSV.write("output/seg$(k)_stats_rcp$(rcp)_noRetreat$(noRetreat).csv",df)
        else
            md = [mode(val[:,i]) for i in 1:size(val)[2]]
            minval = minimum(val,dims=1)'
            maxval = maximum(val,dims=1)'
            df = DataFrame([segNames md minval maxval], :auto)
            names!(df,[:segment,:modalOption,:minOpt,:maxOpt])
            CSV.write("output/segAdapt_stats_rcp$(rcp)_noRetreat$(noRetreat)_$(k).csv",df)
        end
    end
end

function getbricktime(brickfile,tstart,tend)
    years=ncread(brickfile,"time_proj")
    start_ind =findall(x->x==tstart,years)[1]
    end_ind=findall(x-> x==tend,years)[1]
    return start_ind,end_ind
end

# Function: costs as percentage of GDP
function plotMap(tabstr)
    ciamLonLat = CSV.read("data/diva_segment_latlon.csv")
    df = CSV.read(tabstr)


end

# Function: plot costs on map (5-95%, outliers as insets)

# Function: plot dist
# function plotDists(tabstr,globOrSeg="glob")
#     tab = CSV.read(tabstr)
#     if globOrSeg=="glob"
#          ## GMSL Distribution Plot
#          gmslPlot = @df tab density(:gmsl,group=(:brickOutput),label=[:brickOutput])

#     else
#     end

# end
