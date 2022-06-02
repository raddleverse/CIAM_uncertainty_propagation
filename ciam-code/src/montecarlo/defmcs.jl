using Distributions

ciam_default_mcs = @defsim begin
    #slrcost.movefactor  = TriangularDist(0.5,3,1)          # From Diaz (2016); incorporates communication with Mendelsohn and Anthoff and Tol (2014)
    slrcost.movefactor  = Truncated(Normal(1,1),0.5,3)       # From Diaz (2016); incorporates communication with Mendelsohn and Anthoff and Tol (2014)
    slrcost.dvbm        = Truncated(Normal(5.376,2.688),0.0,Inf) # Updated to 2010USD from FUND
    slrcost.vslel       = Truncated(Normal(0.47,0.15),0.0,Inf) # From Viscusi and Aldy (2003)
    slrcost.vslmult     = Truncated(Normal(200,100),0.0,Inf) # From FUND (originally Cline (1992))
    slrcost.wvel        = Truncated(Normal(1.16, 0.46),0.0,Inf) # From Brander et al (2006)
    slrcost.wvpdl       = Truncated(Normal(0.47,0.12),0.0,1.0) # From Brander et al (2006)
end

# hack to use the same interface but only sample a fixed central value
ciam_fixed_mcs = @defsim begin
    slrcost.movefactor  = TriangularDist(1,1,1)          # From Diaz (2016); incorporates communication with Mendelsohn and Anthoff and Tol (2014)
    slrcost.dvbm        = TriangularDist(5.376,5.376,5.376)    # Updated to 2010USD from FUND
    slrcost.vslel       = TriangularDist(0.47,0.47,0.47)     # From Viscusi and Aldy (2003)
    slrcost.vslmult     = TriangularDist(200,200,200)      # From FUND (originally Cline (1992))
    slrcost.wvel        = TriangularDist(1.16,1.16,1.16)      # From Brander et al (2006)
    slrcost.wvpdl       = TriangularDist(0.47,0.47,0.47)     # From Brander et al (2006)
end

function getmcs(vary_ciam)
    if vary_ciam
        return deepcopy(ciam_default_mcs)
    else
        return deepcopy(ciam_fixed_mcs)
    end
end
