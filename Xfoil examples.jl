using Xfoil, Plots, Printf
#read airfoil coordinates from a file
x, y = open("naca2412.dat","r") do f
    x = Float64[]
    y = Float64[]
    for line in eachline(f)
        entries = split(chomp(line))
        push!(x, parse(Float 64, entries[1]))
        push!(y, parse(Float 64, entries[2]))
    end 
    x, y
end
 
#load Xfoil coordinates into Xfoil
Xfoil.set_coordinates(x, y)

#plot the airfoil geometry
scatter(x, y, label="", framestyle=:none, aspect_ratio=10.0, show=true)