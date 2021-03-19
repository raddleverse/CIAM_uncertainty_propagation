## Run all baseline comparisons
# Axes of Variation:
#   RCP, Percentile, NoRetreat, Maintain, Population

# Main Parameters
#brickfile = "/Users/catherineledna/Desktop/BRICK/output_model/BRICK-fastdyn_physical_uniform_26May2020.nc"
using DataFrames # TWmod
brickfile = "/Users/tony/codes/AISFastDyn_Wong_etal_2017/output_model/BRICK-fastdyn_physical_gamma_01Jun2017.nc" # TWmod
ystart=2010
yend=2100
segIDs = false
n=1
choiceDF=DataFrame()

# Load Sea Level Rise
rcp=[85]
pctl = [1,5,50,95,99]
pop = [0,1,2]

# Recreate Delavane Results: Old Population, NoRetreat=F,Fixed=T, GMSL = BRICK
pop=0
noRetreat=false
allowMaintain=false
fixed=true

import Pkg; Pkg.add("Query") # TWmod
Pkg.add("StatsBase") # TWmod
Pkg.add("DataDeps") # TWmod
Pkg.add("Query") # TWmod
using DataDeps # TWmod
using Query # TWmod
using Mimi # TWmod
using NetCDF # TWmod
using StatsBase, CSV # TWmod
include("ciam-code/src/MimiCIAM.jl") # TWmod
include("ciam-code/src/ciamhelper.jl") # TWmod
include("ciam-code/src/brickLSL.jl") # TWmod

#using MimiCIAM # TWmod
m=MimiCIAM.get_model(t=10,fixed=fixed,noRetreat=noRetreat,allowMaintain=allowMaintain)
update_param!(m,:popinput,pop)

# Run  1, 5, 50, 95, 99 - BRICK
(lsl,gmsl,ensInds) = brick_lsl(rcp[1],segIDs,brickfile,n,1,1.01,2010,2100,10,false)
lslr1=lsl[1,:,:]
(lsl,gmsl,ensInds) = brick_lsl(rcp[1],segIDs,brickfile,n,98.99,99.01,2010,2100,10,false)
lslr99=lsl[1,:,:]
(lsl, gmsl,ensInds) = brick_lsl(rcp[1],segIDs,brickfile,n,5,5.01,2010,2100,10,false)
lslr5=lsl[1,:,:]
(lsl, gmsl,ensInds) = brick_lsl(rcp[1],segIDs,brickfile,n,49.95,50.05,2010,2100,10,false)
lslr50=lsl[1,:,:]
(lsl, gmsl,ensInds) = brick_lsl(rcp[1],segIDs,brickfile,n,94.99,95.01,2010,2100,10,false)
lslr95=lsl[1,:,:]
runtag="noRetreatF_maintF_popO"
runAll(lslr1,lslr99,lslr5,lslr50,lslr95,rcp[1],m,runtag)


# Set 2: Flexible
fixed=false
m=MimiCIAM.get_model(t=10,fixed=fixed,noRetreat=noRetreat,allowMaintain=allowMaintain)
runtag="noRetreatF_maintF_popO"
runAll(lslr1,lslr99,lslr5,lslr50,lslr95,rcp[1],m,runtag)


# Set 3: Flexible x PopJones
pop=1
update_param!(m,:popinput,pop)
runtag="noRetreatF_maintF_popJ"
runAll(lslr1,lslr99,lslr5,lslr50,lslr95,rcp[1],m,runtag)

# Set 4: Flexible x PopMerkens
pop=2
update_param!(m,:popinput,pop)
runtag="noRetreatF_maintF_popM"
runAll(lslr1,lslr99,lslr5,lslr50,lslr95,rcp[1],m,runtag)

# Set 5: Flexible x PopJones x NoRetreat
pop=1
update_param!(m,:popinput,pop)
noRetreat=true
update_param!(m,:noRetreat,noRetreat)
runtag="noRetreatT_maintF_popJ"
runAll(lslr1,lslr99,lslr5,lslr50,lslr95,rcp[1],m,runtag)

# Set 6: Flexible x PopJones x AllowMaintain
noRetreat=false
update_param!(m,:noRetreat,noRetreat)
allowMaintain=true
update_param!(m,:allowMaintain,allowMaintain)
runtag="noRetreatF_maintT_popJ"
runAll(lslr1,lslr99,lslr5,lslr50,lslr95,rcp[1],m,runtag)

# Set 7: Flexible x PopJones x NoRetreat x AllowMaintain
noRetreat=true
update_param!(m,:noRetreat,noRetreat)
allowMaintain=true
update_param!(m,:allowMaintain,allowMaintain)
runtag="noRetreatT_maintT_popJ"
runAll(lslr1,lslr99,lslr5,lslr50,lslr95,rcp[1],m,runtag)






function runAll(lsl1,lsl99,lsl5,lsl50,lsl95,rcp,m,runtag)
    update_param!(m,:lslr,lsl1)
    update_param!(m,:rcp,rcp)
    update_param!(m,:percentile,1)
    run(m)
    write_optimal_costs(m;runname=runtag)

    update_param!(m,:lslr,lsl5)
    update_param!(m,:rcp,rcp)
    update_param!(m,:percentile,5)
    run(m)
    write_optimal_costs(m;runname=runtag)

    update_param!(m,:lslr,lsl50)
    update_param!(m,:rcp,rcp)
    update_param!(m,:percentile,50)
    run(m)
    write_optimal_costs(m;runname=runtag)

    update_param!(m,:lslr,lsl95)
    update_param!(m,:rcp,rcp)
    update_param!(m,:percentile,95)
    run(m)
    write_optimal_costs(m;runname=runtag)

    update_param!(m,:lslr,lsl99)
    update_param!(m,:rcp,rcp)
    update_param!(m,:percentile,99)
    run(m)
    write_optimal_costs(m;runname=runtag)
end

function computeOptimalChoicePct(model,year,outdf)

    opt = model[:slrcost,:OptimalOption][year,:]
    lev = model[:slrcost,:OptimalLevel][year,:]
    optlev = opt .* lev
    vals=unique(optlev)
    numSegs=12148

    nameVals=[:Retreat10,:Retreat1,:Retreat100,:Retreat1000,
        :Retreat10000,:NoAdapt,:Protect10,:Protect100,:Protect1000,:Protect10000]

    pcts=[length(findall(x -> x==i,optlev))/numSegs for i in vals]
    temp=DataFrame(pcts')
    names!(temp,nameVals)

    outdf=[outdf;temp]

    return(outdf)

end
