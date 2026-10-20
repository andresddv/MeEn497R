include("Functions.jl")
using .RotorAnalysis
using CCBlade, PyPlot
import .RotorAnalysis: generate_prop_base, run_simulation
# Force a complete reset of any open window to prevent overlapping lines
close("all")
ion()
pygui(true)

# 1. Load airfoil polar data and true 10x7E baseline geometry
af_path = joinpath(dirname(pathof(CCBlade)), "..", "test", "data", "naca4412.dat")
af = AlphaAF("naca4412.dat")
rotor, r, chord, theta = generate_prop_base()

# 2. Set operating simulation conditions (Sweep extended up to J = 0.82)
J_vec = range(0.1, 0.82, length=30) 
Omega = 5400.0 * pi / 30
rho = 1.225

# 3. Compute baseline numerical results for APC 10x7E
eff_base, CT_base, CP_base = run_simulation(rotor, r, chord, theta, af, J_vec, Omega, rho)
CQ_base = CP_base./(2*pi) # Convert CP to CQ for plotting

#=+20% twist modification
theta = theta .* 1.2
eff_mod, CT_mod, CP_mod = run_simulation(rotor, r, chord, theta, af, J_vec, Omega, rho)
CQ_mod = CP_mod./(2*pi) =#

#=chord modification
chord = chord.* 0.8
eff_mod, CT_mod, CP_mod = run_simulation(rotor, r, chord, theta, af, J_vec, Omega, rho)
CQ_mod = CP_mod./(2*pi) =#


#= --- RADIUS SENSITIVITY ANALYSIS ---
# Increase the propeller radius by 10% (e.g., from a 10-inch to an 11-inch propeller)
Rtip_mod = rotor.Rtip * 1.10

# 1. Re-build the rotor structure with the new modified tip radius
rotor_mod = Rotor(rotor.Rhub, Rtip_mod, rotor.B)

# 2. Re-scale the radial stations vector (r) to match the new length
r_mod = (r ./ rotor.Rtip) .* Rtip_mod

# Run simulation with the new modified rotor and scaled radial dimensions
eff_mod, CT_mod, CP_mod = run_simulation(rotor_mod, r_mod, chord, theta, af, J_vec, Omega, rho)
CQ_mod = CP_mod./(2*pi)
=#



# ==============================================================================
# 5. UIUC EXPERIMENTAL DATASET FOR APC 10x7E
# ==============================================================================
exp = [
    0.08  0.111  0.0525;
    0.12  0.109  0.0530;
    0.16  0.106  0.0535;
    0.20  0.103  0.0540;
    0.24  0.099  0.0542;
    0.28  0.095  0.0545;
    0.32  0.090  0.0542;
    0.36  0.085  0.0538;
    0.40  0.079  0.0525;
    0.44  0.072  0.0505;
    0.48  0.065  0.0478;
    0.52  0.057  0.0440;
    0.56  0.048  0.0392;
    0.60  0.040  0.0345;
    0.64  0.031  0.0298;
    0.68  0.023  0.0250;
    0.72  0.015  0.0202;
    0.76  0.008  0.0155;
    0.80  0.002  0.0108;
    0.83  0.000  0.0075
]

J_exp  = exp[:, 1]
CT_exp = exp[:, 2]
CP_exp = exp[:, 3]
CQ_exp = CP_exp./(2*pi) # Convert CP to CQ
eta_exp = (J_exp .* CT_exp) ./ CP_exp


# ==============================================================================
# RADIUS SENSITIVITY PLOTS (COEFFICIENTS & EFFICIENCY)
# ==============================================================================

# Figure 1: Radius Sensitivity Impact on Performance Coefficients (CT, CQ, CP)
figure()
# --- Baseline Configuration (Solid Lines) ---
plot(J_vec, CT_base, "b-", label="CCBlade " * L"C_T" * " (Base\")")
plot(J_vec, CP_base, "r-", label="CCBlade " * L"C_P" * " (Base\")")
plot(J_vec, CQ_base, "g-", label="CCBlade " * L"C_Q" * " (Base\")")

# --- Radius Sensitivity Configuration (Dashed Lines - 11\" Propeller) ---
plot(J_vec, CT_mod, "b--", label="CCBlade " * L"C_T" * " (+10% Radius)")
plot(J_vec, CP_mod, "r--", label="CCBlade " * L"C_P" * " (+10% Radius)")
plot(J_vec, CQ_mod, "g--", label="CCBlade " * L"C_Q" * " (+10% Radius)")

#= --- UIUC Experimental Validation (Markers) ---
plot(J_exp, CT_exp, "bo", label="UIUC Exp " * L"C_T" * " (10x7E)")
plot(J_exp, CP_exp, "ro", label="UIUC Exp " * L"C_P" * " (10x7E)")
plot(J_exp, CQ_exp, "go", label="UIUC Exp " * L"C_Q" * " (10x7E)")=#

xlabel(L"Advance Ratio, $J$")
ylabel("Performance Coefficients")
title("Propeller Geometric Sensitivity: Radius Impact on Coefficients")
grid(true)
legend()
savefig("sensitivity_radius_coefficients.png")


# Figure 2: Radius Sensitivity Impact on Propeller Efficiency (η)
figure()
# --- Efficiency Comparison (Baseline vs Modified vs Experimental) ---
plot(J_vec, eff_base, "g", label="CCBlade " * L"\eta" * " (Base\")")
plot(J_vec, eff_mod, "p--", label="CCBlade " * L"\eta" * " (+10% Radius)")
#=plot(J_exp, eta_exp, "ko", label="UIUC Exp " * L"\eta" * " (10x7E)")=#

xlabel(L"Advance Ratio, $J$")
ylabel(L"Propeller Efficiency, $\eta$")
title("Propeller Geometric Sensitivity: Radius Impact on Efficiency")
grid(true)
legend()
savefig("sensitivity_radius_efficiency.png")

println("Radius sensitivity plots generated successfully!")