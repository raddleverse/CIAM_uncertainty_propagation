##==============================================================================
## processResults.jl
## Downscale BRICK from GMSL to LSL
##
## Original code: Catherine Ledna (18 Mar 2020)
## Modified code: Tony Wong (16 May 2021)
##==============================================================================
using Query
using CSV
using NetCDF
using StatsBase
using DataFrames

##==============================================================================
## Supporting Functions to Downscale BRICK from GMSL to LSL

"""
    get_fingerprints()

Retrieve BRICK fingerprints from NetCDF file - will download the file to a
folder `data` directory
"""
function get_fingerprints()

    fp_dir = joinpath(@__DIR__, "..", "data")
    isdir(fp_dir) || mkpath(fp_dir)
    fp_file = joinpath(fp_dir, "FINGERPRINTS_SLANGEN_Bakker.nc")
    if !isfile(fp_file)
        url = "https://github.com/scrim-network/BRICK/raw/master/fingerprints/FINGERPRINTS_SLANGEN_Bakker.nc"
        download(url, fp_file)
    end

    fplat = ncread(fp_file,"lat")
    fplon = ncread(fp_file,"lon")
    fpAIS = ncread(fp_file,"AIS")
    fpGSIC = ncread(fp_file,"GLAC")
    fpGIS = ncread(fp_file,"GIS")
    ncclose()

    return fplat,fplon,fpAIS,fpGSIC,fpGIS
end

"""
    get_brickGMSL_netcdf(gmslfile::String, rcp::Union{String, Number})

Get brick ensemble members for specified RCP from NetCDF file
Returns time x ens arrays for brick components
"""
function get_brickGMSL_netcdf(gmslfile::String, rcp::Union{String, Number})

    brAIS= ncread(gmslfile,"AIS_RCP$(rcp)")
    brGSIC = ncread(gmslfile,"GSIC_RCP$(rcp)")
    brGIS = ncread(gmslfile,"GIS_RCP$(rcp)")
    brTE = ncread(gmslfile,"TE_RCP$(rcp)")
    brLWS = ncread(gmslfile,"LWS_RCP$(rcp)")
    brGMSL = ncread(gmslfile,"GlobalSeaLevel_RCP$(rcp)")
    btime = ncread(gmslfile,"time_proj")
    ncclose()

    return btime, brAIS, brGSIC, brGIS, brTE, brLWS, brGMSL

end

"""
    get_brickGMSL_rdata(gmslfile::String, rcp::Union{String, Number})

Get brick ensemble members for specified RCP from RData file and return time x ens
arrays for brick components
"""
function get_brickGMSL_rdata(gmslfile::String, rcp::Union{String, Number})

    if !isfile(gmslfile)
        url = "https://zenodo.org/record/3628215/files/sample_projections.RData"
        download(url, gmslfile)
    end

    brick = RData.load(gmslfile)
    brAIS = brick["ais.rcp$(rcp)"]
    brGSIC = brick["gic.rcp$(rcp)"]
    brGIS = brick["gis.rcp$(rcp)"]
    brTE = brick["te.rcp$(rcp)"]
    brGMSL = brick["slr.rcp$(rcp)"]
    brLWS = brGMSL - (brAIS + brGSIC + brGIS + brTE) # LWS not explicitly stored on this file, but it balances the GMSL
    btime = brick["proj.time"]

    return btime, brAIS, brGSIC, brGIS, brTE, brLWS, brGMSL

end

"""
    get_brickGMSL_zip(gmslfile::String, local::rcp::Union{String, Number})

Get brick ensemble members for specified RCP from zip file URL or local path,
and return time x ens arrays for brick components
"""
function get_brickGMSL_zip(gmslfile::String, rcp::Union{String, Number})

    results_dir = joinpath(@__DIR__, "..", "data", "lslr")#TODO HERE NOW!!
    filepath_AIS = joinpath(results_dir,"projections_csv/RCP$(rcp)/projections_antarctic_RCP$(rcp)_sneasybrick_20M_19-02-2022.csv")
    filepath_GSIC = joinpath(results_dir,"projections_csv/RCP$(rcp)/projections_glaciers_RCP$(rcp)_sneasybrick_20M_19-02-2022.csv")
    filepath_GIS = joinpath(results_dir,"projections_csv/RCP$(rcp)/projections_greenland_RCP$(rcp)_sneasybrick_20M_19-02-2022.csv")
    filepath_TE = joinpath(results_dir,"projections_csv/RCP$(rcp)/projections_thermal_RCP$(rcp)_sneasybrick_20M_19-02-2022.csv")
    filepath_GMSL = joinpath(results_dir,"projections_csv/RCP$(rcp)/projections_gmsl_RCP$(rcp)_sneasybrick_20M_19-02-2022.csv")
    filepath_LWS = joinpath(results_dir,"projections_csv/RCP$(rcp)/projections_landwater_storage_sl_RCP$(rcp)_sneasybrick_20M_19-02-2022.csv")
    filepath_MAP = joinpath(results_dir,"projections_csv/RCP$(rcp)/projections_MAP_RCP$(rcp)_sneasybrick_20M_19-02-2022.csv")

    if !isfile(filepath_AIS) | !isfile(filepath_GSIC) | !isfile(filepath_GIS) |
       !isfile(filepath_TE) | !isfile(filepath_GMSL) | !isfile(filepath_LWS) |
       !isfile(filepath_MAP)
        url = "https://zenodo.org/record/6461560/files/sneasybrick_projections_csv.zip"
        download(url, gmslfile)
        run(pipeline(`unzip $(gmslfile) -d $(results_dir)`));
    end

    brAIS = Matrix(CSV.read(filepath_AIS,DataFrame))
    brGSIC = Matrix(CSV.read(filepath_GSIC,DataFrame))
    brGIS = Matrix(CSV.read(filepath_GIS,DataFrame))
    brTE = Matrix(CSV.read(filepath_TE,DataFrame))
    brGMSL = Matrix(CSV.read(filepath_GMSL,DataFrame))
    brLWS = Matrix(CSV.read(filepath_LWS,DataFrame))
    brMAP = CSV.read(filepath_MAP,DataFrame)
    btime = brMAP[!,:YEAR]

    return btime, brAIS, brGSIC, brGIS, brTE, brLWS, brGMSL

end

"""
    get_lonlat(segIDs)

Get CIAM lonlat tuples for specified segIDs, segID order does not matter; will sort
tuples alphabetically by segment name
"""
function get_lonlat(segIDs)

    ciamlonlat = CSV.read(joinpath(@__DIR__,"..","data","diva_segment_latlon.csv"), DataFrame)

    if segIDs == false
        filt = DataFrame(ciamlonlat)
    else
        filt = ciamlonlat |> @filter(_.segid in segIDs) |> Datarame
    end

    sort!(filt, :segments)
    lons = filt.longi
    lats = filt.lati

    return collect(zip(lons,lats))
end

"""
    choose_ensemble_members(time, ens, n, low, high, yend, ensInds)

Choose n ensemble members randomly within specified percentile range
if percentile range is one number, will return that percentile (with respect to GMSL in specified year)
time - BRICK time vector from netcdf
ens - BRICK GMSL matrix, time x num ensembles
low - minimum percentile threshold (integer - e.g. 5 = 5th percentile)
high - maximum percentile threshold
"""
function choose_ensemble_members(time, ens, n, low, high, yend, ensInds)

    if length(time)==size(ens)[1]
        end_year = findall(x -> x==yend, time)[1]

        # want to ignore NaNs later
        idx_good = findall(!isnan, ens[end_year,:])

        if low==high
            # pick just the one ensemble member that fits,
            # but get rid of NANs first
            val_low = percentile(ens[end_year,idx_good],low) # only need one because we know they're equal here
            # then, see if there's a single member of the ensemble that has the desired percentile
            idx_perc = findall(x -> x==val_low, ens[end_year,idx_good]) # idx_perc is idx within ens_inds
            if length(idx_perc) > 1
                idx_perc=idx_perc[1] # if more than 1, just take the first one
            elseif length(idx_perc) == 0
                # take the nearest one
                idx_perc = findmin(broadcast(abs, ens[end_year, idx_good] .- val_low))[2]
                return idx_good[idx_perc]
            else
                return idx_good[idx_perc]
            end
        else
            if ensInds==false
                # if you want to restrict to only between certain quantiles within the ensemble
                val_low = percentile(ens[end_year,idx_good],low)
                val_high = percentile(ens[end_year,idx_good],high)
                # pick ensemble members
                ens_inds = findall(x -> x >= val_low && x <= val_high, ens[end_year,:])
                if length(ens_inds)>1
                    chosen_inds = ens_inds[sample(1:end,n,replace=false)]
                else
                    chosen_inds= ens_inds
                end
                return chosen_inds
            else
                chosen_inds = ensInds
                return chosen_inds
            end
        end

    else
        println("Error: time dimension mismatch")
    end
end

"""
    downscale_brick(brickcomps,lonlat, ensInds, ystart=2010, yend=2100, tstep=10)
Downscale BRICK gmsl to lsl for all segments and ensembles of interest with the
input arguments:

brickcomps - BRICK components (time x ens matrices corresponding to brick gmsl components)
lonlat - vector of (lon,lat) tuples, sorted corresp to segment name alphabetical order
Output:
lsl_out: ens x time x segment array of local sea levels, sorted in alphabetical order by segment name
GMSL: global mean sea levels corresponding to local sea level arrays (time x ens)
"""
function downscale_brick(brickcomps, lonlat, ensInds, ystart=2010, yend=2100, tstep=10)
    # To do - check with vectors of lat, lon
    (fplat,fplon,fpAIS,fpGSIC,fpGIS) = get_fingerprints()
    (btime,AIS,GSIC,GIS,TE,LWS,GMSL) = brickcomps

    # Select indices of time of interest, with respect to timestep
    tinds = findall(x -> x .>= ystart && x .<=yend, btime)
    years = collect(ystart:yend)
    yinds = findall(x -> x % tstep==0, years)
    # Need to normalize LSL relative to 2000
    inorm = findall(x -> x==2000, btime)

    tdim=length(btime)

    if length(years)==length(tinds)
        tinds = tinds[yinds]
    else
        println("Error: years outside of bounds")
        return nothing
    end

    num_ens = length(ensInds)

    # Output matrix: ens x time x segment
    lsl_out = zeros(num_ens, length(tinds), length(lonlat))

    # Trim component vectors to timesteps and ensembles. Assume interval is 1 year
    if tdim==size(AIS)[1] # check that time dimension is 1
        # for normalizing
        AIS_norm = AIS[inorm,ensInds]
        GSIC_norm = GSIC[inorm,ensInds]
        GIS_norm = GIS[inorm,ensInds]
        TE_norm = TE[inorm,ensInds]
        LWS_norm = LWS[inorm,ensInds]
        GMSL_norm = GMSL[inorm,ensInds]
        # actual projections
        AIS = AIS[tinds,ensInds]
        GSIC = GSIC[tinds,ensInds]
        GIS = GIS[tinds,ensInds]
        TE = TE[tinds,ensInds]
        LWS = LWS[tinds,ensInds]
        GMSL = GMSL[tinds,ensInds]
    else
        println("Error: time dimension is not 1 for brick components")
        return nothing
    end

    for f in 1:length(lonlat) # Loop through lonlat tuples

        lon = lonlat[f][1]
        lat = lonlat[f][2]
        # Convert Longitude to degrees East
        # CIAM Lat is already in (-90,90) by default
        if lon <0
            lon = lon + 360
        end

        # Find fingerprint degrees nearest to lat,lon
        ilat = findall(isequal(minimum(abs.(fplat.-lat))),abs.(fplat.-lat))
        ilon = findall(isequal(minimum(abs.(fplon.-lon))),abs.(fplon.-lon))


        # Take average of closest lat/lon values
        fpAIS_flat = collect(skipmissing(Iterators.flatten(fpAIS[ilon,ilat])))
        fpGSIC_flat = collect(skipmissing(Iterators.flatten(fpGSIC[ilon,ilat])))
        fpGIS_flat = collect(skipmissing(Iterators.flatten(fpGIS[ilon,ilat])))

        fpAIS_loc = mean(fpAIS_flat[isnan.(fpAIS_flat).==false],dims=1)[1]
        fpGSIC_loc = mean(fpGSIC_flat[isnan.(fpGSIC_flat).==false],dims=1)[1]
        fpGIS_loc = mean(fpGIS_flat[isnan.(fpGIS_flat).==false],dims=1)[1]
        fpTE_loc = 1.0
        fpLWS_loc=1.0

        # Keep searching nearby lat/lon values if fingerprint value is NaN unless limit is hit
        inc =1

        while isnan(fpAIS_loc) || isnan(fpGIS_loc) || isnan(fpGSIC_loc) && inc<5

            newlonStart = lon_subtractor.(fplon[ilon],inc)[1]
            newlatStart = lat_subtractor.(fplat[ilat],inc)[1]
            newlonEnd = lon_adder.(fplon[ilon],inc)[1]
            newlatEnd = lat_adder.(fplat[ilat],inc)[1]

            latInd1 = minimum(findall(isequal(minimum(abs.(fplat.-newlatStart))),abs.(fplat.-newlatStart)))
            #minimum(findall(x-> x in newlatStart,fplat))
            latInd2 = maximum(findall(isequal(minimum(abs.(fplat.-newlatEnd))),abs.(fplat.-newlatEnd)))
            #maximum(findall(x -> x in newlatEnd,fplat))

            lonInd1 = minimum(findall(isequal(minimum(abs.(fplon.-newlonStart))),abs.(fplon.-newlonStart)))
            #minimum(findall(x-> x in newlonStart,fplon))
            lonInd2 = maximum(findall(isequal(minimum(abs.(fplon.-newlonEnd))),abs.(fplon.-newlonEnd)))
            #maximum(findall(x -> x in newlonEnd,fplon))

            if latInd2 < latInd1
                latInds=[latInd1; 1:latInd2]
            else
                latInds=latInd1:latInd2
            end

            if lonInd2 < lonInd1
                lonInds=[lonInd1; 1:lonInd2]
            else
                lonInds = lonInd1:lonInd2
            end

            fpAIS_flat = collect(skipmissing(Iterators.flatten(fpAIS[lonInds,latInds])))
            fpGSIC_flat = collect(skipmissing(Iterators.flatten(fpGSIC[lonInds,latInds])))
            fpGIS_flat = collect(skipmissing(Iterators.flatten(fpGIS[lonInds,latInds])))

            fpAIS_loc = mean(fpAIS_flat[isnan.(fpAIS_flat).==false],dims=1)[1]
            fpGSIC_loc = mean(fpGSIC_flat[isnan.(fpGSIC_flat).==false],dims=1)[1]
            fpGIS_loc = mean(fpGIS_flat[isnan.(fpGIS_flat).==false],dims=1)[1]

            inc = inc + 1

        end

        # If still NaN, throw an error
        if isnan(fpAIS_loc) || isnan(fpGIS_loc) || isnan(fpGSIC_loc)
            println("Error: no fingerprints found for ($(lon),$(lat))")
            return nothing
        end

       # Multiply fingerprints by BRICK ensemble members
       if ndims(AIS) > 1
            for n in 1:size(AIS)[2] # loop through ensemble members
                lsl_out[n, :, f] = fpGIS_loc * GIS[:,n] + fpAIS_loc * AIS[:,n] + fpGSIC_loc * GSIC[:,n] +
                                   fpTE_loc * TE[:,n] + fpLWS_loc * LWS[:,n]
                # CIAM - LSL should be sea-level change relative to year 2000
                lsl_norm = fpGIS_loc * GIS_norm[n] + fpAIS_loc * AIS_norm[n] + fpGSIC_loc * GSIC_norm[n] +
                           fpTE_loc * TE_norm[n] + fpLWS_loc * LWS_norm[n]
                lsl_out[n, :, f] = lsl_out[n, :, f] .- lsl_norm
            end
        else
            lsl_out[1, :, f] = fpGIS_loc * GIS[:] + fpAIS_loc * AIS[:] + fpGSIC_loc * GSIC[:] +
                fpTE_loc * TE[:] + fpLWS_loc * LWS[:]
            # CIAM - LSL should be sea-level change relative to year 2000
            lsl_norm = fpGIS_loc * GIS_norm + fpAIS_loc * AIS_norm + fpGSIC_loc * GSIC_norm +
                       fpTE_loc * TE_norm + fpLWS_loc * LWS_norm
            lsl_out[1, :, f] = lsl_out[1, :, f] .- lsl_norm
        end

    end # End lonlat tuple

    return lsl_out,GMSL
end

"""
    brick_lsl(rcp,segIDs,brickfile,n,low=5,high=95,ystart=2010,yend=2100,tstep=10,ensInds=false)

Driver function to downscale BRICK gmsl for specified segments
"""
function brick_lsl(rcp,segIDs,brickfile,n,low=5,high=95,ystart=2010,yend=2100,tstep=10,ensInds=false)
    # HERE - if you want to use a different set of SLR projections, or projections
    # that are stored in a different format, a new get_brickGMSL_xxx might be needed
    #brickGMSL = get_brickGMSL_rdata(brickfile,rcp)
    brickGMSL = get_brickGMSL_zip(brickfile, rcp)
    brickEnsInds = choose_ensemble_members(brickGMSL[1],brickGMSL[7],n,low,high,yend,ensInds)
    lonlat = get_lonlat(segIDs)
    (lsl,gmsl) = downscale_brick(brickGMSL, lonlat, brickEnsInds, ystart, yend, tstep)

    return lsl,gmsl,brickEnsInds
end

##==============================================================================
## Small Helper Functions

function adder(maxval)
    function y(point,n)
        if point + n > maxval
            return point + n - maxval
        else
            return point + n
        end
    end
end

function subtractor(minval,maxval)
    function y(point,n)
        if point - n < minval
            return min(maxval,point - n + maxval)
        else
            return point - n
        end
    end
end

lon_subtractor = subtractor(1,360)
lon_adder = adder(360)
lat_adder = adder(180)
lat_subtractor = subtractor(1,180)
