using Xfoil, Plots, Printf
#read airfoil coordinates from a file
x, y = open("naca2412.dat","r") do f
    x = Float64[]
    y = Float64[]
    for line in eachline(f)
        entries = split(chomp(line))
        push!(x, parse(Float64, entries[1]))
        push!(y, parse(Float64, entries[2]))
    end 
    x, y
end
 
#load Xfoil coordinates into Xfoil
Xfoil.set_coordinates(x, y)

#plot the airfoil geometry
f = scatter(x, y, label="", framestyle=:none, aspect_ratio=1.0, show=true)
display(f)

# repanel using XFOIL's `PANE` command
xr, yr = Xfoil.pane()

# plot the refined airfoil geometry
j = scatter(xr, yr, label="", framestyle=:none, aspect_ratio=1.0, show=true)
display(j)

# set operating conditions
alpha = -9:1:14 # range of angle of attacks, in degrees
re = 1e5 # Reynolds number

# initialize outputs
n_a = length(alpha)
c_l = zeros(n_a)
c_d = zeros(n_a)
c_dp = zeros(n_a)
c_m = zeros(n_a)
converged = zeros(Bool, n_a)

# determine airfoil coefficients across a range of angle of attacks
for i = 1:n_a
    c_l[i], c_d[i], c_dp[i], c_m[i], converged[i] = Xfoil.solve_alpha(alpha[i], re; iter=100, reinit=true)
end

# print results
println("Angle\t\tCl\t\tCd\t\tCm\t\tConverged")
for i = 1:n_a
    @printf("%8f\t%8f\t%8f\t%8f\t%d\n",alpha[i],c_l[i],c_d[i],c_m[i],converged[i])
end

# plot results
p1 = plot(alpha, c_l,
    label="",
    xlabel="Angle of Attack (degrees)",
    ylabel="Lift Coefficient")

p2 = plot(alpha, c_d,
    label="",
    xlabel="Angle of Attack (degrees)",
    ylabel="Drag Coefficient")

p3 = plot(alpha, c_m,
    label="",
    xlabel="Angle of Attack (degrees)",
    ylabel="Moment Coefficient")

display(p1)

display(p2)

display(p3)