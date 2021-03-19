### Make Figures 
using StatsPlots
using DataFrames

function run_fastdynamics()

end


function plotDists(table, cols)
    

end

### Updated Plots
w3 = CSV.read("../test/FDcombined.csv")
noFD = CSV.read("test/noFD_wRetreat.csv")
noFDnr = CSV.read("test/noFD_noRetreat.csv")
names!(noFDnr,[:NPV,:WetlandLoss,:DrylandLoss,:StormLoss,:GMSL])

noFD = stack(noFD)
noFDnr = stack(noFDnr)
noFDnr[:type]="NoFD-NoRetreat"
gmslNoFD = @from i in noFD begin
    @where i.variable==:GMSL
    @select {i.variable,i.value,i.type} 
    @collect DataFrame
end
gmslNoFDv = collect(gmslNoFD[:value])
gmslFD = collect(gmsl[:value])

npv = @from i in w3 begin
    @where i.variable==:NPV 
    @select {i.variable,i.value,i.type} 
    @collect DataFrame
end
npv[:gmsl] = gmslFD

npvNoFD = @from i in noFD begin
    @where i.variable==:NPV
    @select {i.variable,i.value,i.type} 
    @collect DataFrame
end
npvNoFD[:gmsl] = gmslNoFDv
npvAll = [npv;npvNoFD]

npvNoFDNr = @from i in noFDnr begin
    @where i.variable==:NPV
    @select {i.variable,i.value,i.type} 
    @collect DataFrame
end
npvNoFDNr[:gmsl] = gmslNoFDv
npvAll = [npvAll;npvNoFDNr]

wetland = @from i in w3 begin
    @where i.variable==:WetlandLoss
    @select {i.variable,i.value,i.type} 
    @collect DataFrame
end

dryland = @from i in w3 begin
    @where i.variable==:DrylandLoss 
    @select {i.variable,i.value,i.type} 
    @collect DataFrame
end
dryland2 = @from i in noFD begin
    @where i.variable==:Dryland 
    @select {i.variable,i.value,i.type} 
    @collect DataFrame
end
drylandAll = [dryland;dryland2]

storm = @from i in w3 begin
    @where i.variable==:StormLoss
    @select {i.variable,i.value,i.type} 
    @collect DataFrame
end

gmsl = @from i in w3 begin
    @where i.variable==:GMSL 
    @select {i.variable,i.value,i.type} 
    @collect DataFrame
end

gmsl2 = @from i in w3 begin
    @where i.variable==:GMSL && i.type=="FD-Retreat"
    @select {i.variable,i.value,i.type} 
    @collect DataFrame
end


plot(
   l1,l2,l3,layout=(3,1)
)
l1 = @df npvAll density(:value,group=(:type),label=["Fast Dynamics - No Retreat","Fast Dynamics - With Retreat","No Fast Dynamics - With Retreat","No Fast Dynamics - No Retreat"])
xlabel!(l1,"NPV (Billion 2010USD)")
savefig(l1,"test/plots/NPVwFDcomp3.png")

wet2 = @from i in wetland begin
    @where i.type=="FD-Retreat"
    @select {i.variable,i.value,i.type} 
    @collect DataFrame
end

l2 = @df wet2 density(:value)
xlabel!(l2,"Wetland Loss (Km^2)")

gmslAll = [gmsl2;gmslNoFD]
l4 = @df gmslAll density(:value,group=(:type),label=["Fast Dynamics","No Fast Dynamics"])
xlabel!(l4,"GMSL, 2100 (m)")
savefig(l4,"test/plots/GMSLcomp.png")


l3 = @df drylandAll density(:value,group=(:type),label=["Fast Dynamics - No Retreat","Fast Dynamics - With Retreat","No Fast Dynamics - With Retreat"])
xlabel!(l3, "Dryland Loss, (Km^2)")
savefig(l3,"test/plots/DrylandWFDcomp2.png")

l4 = @df storm density(:value,group=(:type),label=["Fast Dynamics - No Retreat","Fast Dynamics - With Retreat"])
xlabel!(l4,"Expected Storm Mortality")
l5 = @df gmsl density(:value,group=(:variable))
xlabel!(l5, "GMSL 2100 (m)")

# Covariate plot 
npvWR = @from i in npv begin
    @where i.type==:"FD-Retreat"
    @select {i.value} 
    @collect DataFrame
end
npvWR[:gmsl] = gmsl[:valGMSL]


npvNR = @from i in npv begin
    @where i.type==:"FD-NoRetreat"
    @select {i.value} 
    @collect DataFrame
end
npvNR[:gmsl] = gmsl[:valGMSL]


c1 = @df npvWR corrplot(cols(1:2), grid = false)
M = [collect(gmsl.valGMSL) fd_wRetreat[1]]
c1 = cornerplot(M, label=["GMSL 2100 (m) - Fast Dynamics","NPV - With Retreat"],compact=true)
savefig(c1,"test/plots/NPVwFDcorrplot.png")


c2 = @df npvNR corrplot(cols(1:2), grid = false)
M2 = [collect(gmsl.valGMSL) fd_noRetreat[1]]
c2 = cornerplot(M2, label=["GMSL 2100 (m) - Fast Dynamics","NPV - No Retreat"],compact=true)
savefig(c2,"test/plots/NPVwFDNrcorrplot.png")

drylandR_WFD = @from i in dryland begin
    @where i.type=="FD-Retreat" 
    @select {i.variable,i.value,i.type} 
    @collect DataFrame
end
drylandNR_WFD = @from i in dryland begin
    @where i.type=="FD-NoRetreat" 
    @select {i.variable,i.value,i.type} 
    @collect DataFrame
end


gmslAll2 = [gmsl;gmslNoFD]
M3 = [collect(gmsl2.value) collect(drylandNR_WFD.value)]
c3 = cornerplot(M3,label=["GMSL 2100 (m) - Fast Dynamics","Dryland Loss (Km^2) - No Retreat"],compact=true)




M2 = [y[5][10,:] y[1]]
c2 = cornerplot(M2, label=["GMSL 2100","NPV"], compact=true)

corrplots=plot(
    c1,c2,layout=(1,2)
)
savefig(corrplots,joinpath(@__DIR__,"../plots/NPV_vs_GMSL.png","width",10,"height",5)

lo = percentile(seg_npv,5)
hi = percentile(seg_npv,95)
slr_inds = findall(x -> x>0, seg_lsl)
inds = findall(x-> x<=hi && x >=lo, seg_npv)
fin_inds = intersect(slr_inds,inds)

lsl_filt = seg_lsl[fin_inds]
seg_npv_filt=seg_npv[fin_inds]
M3 = [lsl_filt seg_npv_filt]

c3 = cornerplot(M3,label=["LSL 2100","Segment NPV"],compact=true,markercolor="blue")