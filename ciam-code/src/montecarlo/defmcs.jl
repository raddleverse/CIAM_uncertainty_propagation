using Distributions 

ciam_default_mcs = @defsim begin
    movefactor= TriangularDist(0.5,3,1)          # From Diaz (2016); incorporates communication with Mendelsohn and Anthoff and Tol (2014)
    dvbm= Truncated(Normal(5.376,2.688),0.0,Inf) # Updated to 2010USD from FUND 
    vslel = Truncated(Normal(0.47,0.15),0.0,Inf) # From Viscusi and Aldy (2003)
    vslmult = Truncated(Normal(200,100),0.0,Inf) # From FUND (originally Cline (1992))
    wvel = Truncated(Normal(1.16, 0.46),0.0,Inf) # From Brander et al (2006)
    wvpdl = Truncated(Normal(0.47,0.12),0.0,1.0) # From Brander et al (2006)
end

function getmcs()
    return deepcopy(ciam_default_mcs)
end