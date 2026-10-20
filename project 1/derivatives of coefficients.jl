using Xfoil, Printf

# read airfoil into XFOIL
open("naca2412.dat", "r") do f
    x = Float64[]
    y = Float64[]
    for line in eachline(f)
        entries = split(chomp(line))
        push!(x, parse(Float64, entries[1]))
        push!(y, parse(Float64, entries[2]))
    end
    Xfoil.set_coordinates(x,y)
end

# repanel using XFOIL's `PANE` command
Xfoil.pane()

# set operating conditions
alpha = -9:1:14
re = 1e5
mach = 0.0

# set step size
h = 1e-6

# initialize outputs
n_a = length(alpha)
c_l_a = zeros(n_a)
c_d_a = zeros(n_a)
c_dp_a = zeros(n_a)
c_m_a = zeros(n_a)
converged = zeros(Bool, n_a)

for i = 1:n_a
    c_l1, c_d1, c_dp1, c_m1, converged[i] = Xfoil.solve_alpha(alpha[i], re; mach, iter=100, reinit=true)
    c_l2, c_d2, c_dp2, c_m2, converged[i] = Xfoil.solve_alpha(alpha[i]+h, re; mach, iter=100, reinit=true)
    c_l_a[i] = (c_l2 - c_l1)/h * 180/pi
    c_d_a[i] = (c_d2 - c_d1)/h * 180/pi
    c_m_a[i] = (c_m2 - c_m1)/h * 180/pi
end

# print results
println("Angle\t\tdClda\t\tdCdda\t\tdCmda\t\tConverged")
for i = 1:n_a
  @printf("%8f\t%8f\t%8f\t%8f\t%d\n",alpha[i],c_l_a[i],c_d_a[i],c_m_a[i],converged[i])
end