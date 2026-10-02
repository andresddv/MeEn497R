

#######
# Compare lift-to-drag curves for several thicknesses.
function run_xfoil(x0, y0, alpha, re, thickness_factor)
    x_new = copy(x0)
    y_new = copy(y0)
    for i in 1:fld(length(y_new), 2)
        j = length(y_new) - i + 1
        y_bar = (y_new[i] + y_new[j]) / 2
        t = (y_new[i] - y_new[j]) * thickness_factor
        y_new[i] = y_bar + t / 2
        y_new[j] = y_bar - t / 2
    end
    Xfoil.set_coordinates(x_new, y_new)
    n = length(alpha)
    cl = zeros(n)
    cd = zeros(n)
    cdp = zeros(n)
    cm = zeros(n)
    conv = zeros(Bool, n)
    for i in eachindex(alpha)
        cl[i], cd[i], cdp[i], cm[i], conv[i] =
            Xfoil.solve_alpha(alpha[i], re; iter=100, reinit=true)
    end
    return cl, cd, cdp, cm, conv
end

p5 = plot(xlabel="Lift coefficient", ylabel="Lift-to-drag ratio", xlims=(-10, 10), ylims=(-500, 500))
for (factor, color) in [(0.7, :black), (0.9, :green), (1.1, :pink),
                        (1.2, :orange), (1.5, :red)]
    cl_t, cd_t, _, _, _ = run_xfoil(x, y, alpha, re, factor)
    plot!(p5, cl_t, cl_t ./ cd_t,
        label="t = $(round(Int, 10 * factor))%", color=color)
end
display(p5)


c_t = [(0.8, :black, "t = 8%"), 
        (0.9, :green, "t = 9%"), (1.1, :pink, "t = 11%"), (1.2, :orange, "t = 12%"), (1.5, :red, "t = 15%")]

for (i, clr, lb) in c_t
    c_l, c_d, c_dp, c_m = run_xfoil(x, y, re, i)
    plot!(p5, lift_over_drag, alpha,
        label=lb, color=clr)
end

display(p5)