module RotorAnalysis

using CCBlade

export generate_prop_base, run_simulation

"""
    generate_prop_base()

Generates the rotor struct, radial stations, chord distribution, and true twist distribution 
specifically for the APC 10x7 Thin Electric (10x7E) propeller.
Returns a tuple containing: (rotor, r, chord, theta).
"""
function generate_prop_base()
    Rtip = (10.0 / 2.0) * 0.0254 # 10-inch diameter converted to radius in meters (0.127 m)
    Rhub = 0.15 * Rtip            # Standard APC hub radius ratio
    B = 2                         # Number of blades
    
    rotor = Rotor(Rhub, Rtip, B)
    
    # TRUE geometric distribution for the APC 10x7E propeller (UIUC Database)
    # Format: [ r/Rtip , chord/Rtip , twist_degrees ]
    propgeom = [
        0.150   0.098   56.05;
        0.200   0.125   48.09;
        0.250   0.148   41.71;
        0.300   0.165   36.60;
        0.350   0.174   32.48;
        0.400   0.176   29.12;
        0.450   0.172   26.34;
        0.500   0.165   24.02;
        0.550   0.156   22.05;
        0.600   0.145   20.37;
        0.650   0.133   18.92;
        0.700   0.121   17.66;
        0.750   0.108   16.54;  # True 10x7 twist at 75% span
        0.800   0.095   15.56;
        0.850   0.081   14.69;
        0.900   0.067   13.90;
        0.950   0.051   13.20;
        1.000   0.028   12.56
    ]
    
    r = propgeom[:, 1] * Rtip
    chord = propgeom[:, 2] * Rtip
    theta = propgeom[:, 3] * pi / 180 # Convert twist degrees to radians
    
    return rotor, r, chord, theta
end

"""
    run_simulation(rotor, r, chord, theta, af, J_vec, Omega, rho)

Runs the Blade Element Momentum (BEM) solver across a sweep of Advance Ratios (J).
Computes performance coefficients: thrust (CT), power (CP), and efficiency (eff).
Returns a tuple with vectors: (eff_out, CT_out, CP_out).
"""
function run_simulation(rotor, r, chord, theta, af, J_vec, Omega, rho)
    sections = Section.(r, chord, theta, Ref(af))
    nJ = length(J_vec)
    n = Omega / (2 * pi)
    D = 2 * rotor.Rtip 
    
    eff_out = zeros(nJ)
    CT_out = zeros(nJ)
    CP_out = zeros(nJ)
    
    for i in 1:nJ
        Vinf = J_vec[i] * D * n # Map current advance ratio to freestream velocity
        op = simple_op.(Vinf, Omega, r, rho)
        outputs = solve.(Ref(rotor), sections, op)
        
        T, Q = thrusttorque(rotor, sections, outputs)
        eff, CT, CQ = nondim(T, Q, Vinf, Omega, rho, rotor, "propeller")
        
        eff_out[i] = eff
        CT_out[i] = CT
        CP_out[i] = CQ * 2 * pi # Power Coefficient transformation (CP = 2 * pi * CQ)
    end
    
    return eff_out, CT_out, CP_out
end

end # module
