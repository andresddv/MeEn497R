using Xfoil, Plots, Printf

# 1. Leer coordenadas del perfil (NACA 1410)
x, y = open("naca1410.dat","r") do f
    x_coords = Float64[]
    y_coords = Float64[]
    for line in eachline(f)
        entries = split(chomp(line))
        if length(entries) >= 2 && tryparse(Float64, entries[1]) !== nothing
            push!(x_coords, parse(Float64, entries[1]))
            push!(y_coords, parse(Float64, entries[2]))
        end
    end 
    x_coords, y_coords
end
 
Xfoil.set_coordinates(x, y)
xr, yr = Xfoil.pane()

# 2. Configuración de simulación XFOIL interna
alpha = -16:0.25:16 

# Inicializar los 3 gráficos independientes e individuales con leyendas transparentes
p1 = plot(xlabel="Angle of Attack (degrees)", ylabel="Lift Coefficient", legend_background_color=:transparent, legend_foreground_color=:transparent)
p2 = plot(xlabel="Angle of Attack (degrees)", ylabel="Drag Coefficient", legend_background_color=:transparent, legend_foreground_color=:transparent)
p3 = plot(xlabel="Angle of Attack (degrees)", ylabel="Moment Coefficient", legend_background_color=:transparent, legend_foreground_color=:transparent)

nasa_data = [(3e6, :blue, "3e6"), 
              (6e6, :orange, "6e6")]

# Data extracted from NASA report
nasa_3m_alpha = [-14.0, -12.0, -10.0, -8.0, -6.0, -4.0, -2.0,  0.0,  2.0,  4.0,  6.0,  8.0, 10.0, 12.0, 13.5]
nasa_3m_cl    = [-0.83, -0.74, -0.61, -0.44, -0.23, -0.02,  0.19,  0.39,  0.59,  0.79,  0.98,  1.16,  1.30,  1.38,  1.14]
nasa_3m_cd    = [0.0175, 0.0125, 0.0098, 0.0078, 0.0066, 0.0062, 0.0061, 0.0062, 0.0066, 0.0075, 0.0088, 0.0110, 0.0142, 0.0195, 0.0270]
nasa_3m_cm    = [-0.010, -0.011, -0.011, -0.011, -0.012, -0.012, -0.012, -0.013, -0.013, -0.013, -0.014, -0.014, -0.015, -0.016, -0.024]

nasa_6m_alpha = [-14.0, -12.0, -10.0, -8.0, -6.0, -4.0, -2.0,  0.0,  2.0,  4.0,  6.0,  8.0, 10.0, 12.0, 14.0]
nasa_6m_cl    = [-0.85, -0.76, -0.62, -0.44, -0.23, -0.02,  0.19,  0.39,  0.59,  0.79,  0.98,  1.17,  1.32,  1.42,  1.18]
nasa_6m_cd    = [0.0160, 0.0115, 0.0090, 0.0072, 0.0061, 0.0058, 0.0057, 0.0058, 0.0062, 0.0070, 0.0083, 0.0102, 0.0130, 0.0175, 0.0250]
nasa_6m_cm    = [-0.010, -0.011, -0.011, -0.011, -0.012, -0.012, -0.012, -0.013, -0.013, -0.013, -0.014, -0.014, -0.015, -0.016, -0.025]

for (re, clr, lb) in nasa_data
    # 1. EJECUTAR XFOIL DINÁMICO PARA EL REYNOLDS ACTUAL
    results = [Xfoil.solve_alpha(a, re; iter=150, reinit=true) for a in alpha]
    local converged = [r[5] for r in results]
    
    # Inicializar vectores limpios con la convergencia base de XFOIL
    c_l_clean = [converged[i] ? results[i][1] : NaN for i in 1:length(alpha)]
    c_d_clean = [converged[i] ? results[i][2] : NaN for i in 1:length(alpha)]
    c_m_clean = [converged[i] ? results[i][4] : NaN for i in 1:length(alpha)]
    #=
    # Filtro estricto de estabilidad para limpiar picos de ruido matemático
    for i in 2:length(alpha)
        if !isnan(c_l_clean[i]) && !isnan(c_l_clean[i-1]) && abs(c_l_clean[i] - c_l_clean[i-1]) > 0.15;  c_l_clean[i] = NaN; end
        if !isnan(c_d_clean[i]) && !isnan(c_d_clean[i-1]) && (c_d_clean[i] > 0.045 || abs(c_d_clean[i] - c_d_clean[i-1]) > 0.015); c_d_clean[i] = NaN; end
        if !isnan(c_m_clean[i]) && !isnan(c_m_clean[i-1]) && abs(c_m_clean[i] - c_m_clean[i-1]) > 0.015;  c_m_clean[i] = NaN; end
    end
    =#
    # 2. GRAFICAR LAS PREDICCIONES DE XFOIL (Líneas sólidas continuas)
    plot!(p1, alpha, c_l_clean, label="XFOIL ($lb)", color=clr, linewidth=1.5, legend_title="Source (Re)")
    plot!(p2, alpha, c_d_clean, label="XFOIL ($lb)", color=clr, linewidth=1.5, legend_title="Source (Re)")
    plot!(p3, alpha, c_m_clean, label="XFOIL ($lb)", color=clr, linewidth=1.5, legend_title="Source (Re)")

    # 3. GRAFICAR LOS DATOS DE LA NASA RESPECTIVOS (Líneas punteadas)
    if re == 3e6
        plot!(p1, nasa_3m_alpha, nasa_3m_cl, label="NASA ($lb)", color=clr, linestyle=:dash, linewidth=1.5)
        plot!(p2, nasa_3m_alpha, nasa_3m_cd, label="NASA ($lb)", color=clr, linestyle=:dash, linewidth=1.5)
        plot!(p3, nasa_3m_alpha, nasa_3m_cm, label="NASA ($lb)", color=clr, linestyle=:dash, linewidth=1.5)
    elseif re == 6e6
        plot!(p1, nasa_6m_alpha, nasa_6m_cl, label="NASA ($lb)", color=clr, linestyle=:dash, linewidth=1.5)
        plot!(p2, nasa_6m_alpha, nasa_6m_cd, label="NASA ($lb)", color=clr, linestyle=:dash, linewidth=1.5)
        plot!(p3, nasa_6m_alpha, nasa_6m_cm, label="NASA ($lb)", color=clr, linestyle=:dash, linewidth=1.5)
    end
end

#=
# FUNCIÓN AUXILIAR CORREGIDA: Lee y LIMPIA el ruido de los archivos de Airfoil Tools
function leer_datos_archivo_limpios(nombre_archivo)
    a_file, cl_file, cd_file, cm_file = Float64[], Float64[], Float64[], Float64[]
    if isfile(nombre_archivo)
        open(nombre_archivo, "r") do f
            for line in eachline(f)
                entries = split(chomp(line))
                if length(entries) == 7 && tryparse(Float64, entries[1]) !== nothing
                    push!(a_file,  parse(Float64, entries[1])) # alpha
                    push!(cl_file, parse(Float64, entries[2])) # CL
                    push!(cd_file, parse(Float64, entries[3])) # CD
                    push!(cm_file, parse(Float64, entries[5])) # CM
                end
            end
        end
        
        # APLICAR FILTRO ANTIRRUIDO DINÁMICO A LOS DATOS LEÍDOS
        for i in 2:length(a_file)
            # Filtro para Lift (Cl) de Airfoil Tools
            if abs(cl_file[i] - cl_file[i-1]) > 0.15
                cl_file[i] = NaN
            end
            
            # Filtro para Drag (Cd) de Airfoil Tools (Umbral máximo y saltos)
            if cd_file[i] > 0.045 || abs(cd_file[i] - cd_file[i-1]) > 0.015
                cd_file[i] = NaN
            end
            
            # Filtro para Moment (Cm) de Airfoil Tools
            if abs(cm_file[i] - cm_file[i-1]) > 0.015
                cm_file[i] = NaN
            end
        end
    end
    return a_file, cl_file, cd_file, cm_file
end

# 3. BUCLE PRINCIPAL UNIFICADO
air_tools = [
    (5e4, :black,  "50,000",    "re5e4.txt"), 
    (1e5, :green,  "100,000",   "re1e5.txt"), 
    (2e5, :pink,   "200,000",   "re2e5.txt"), 
    (5e5, :orange, "500,000",   "re5e5.txt"), 
    (1e6, :red,    "1,000,000", "re1e6.txt")
]

for (re, clr, lb, nombre_txt) in air_tools
    # --- PARTE A: Simulación interactiva con XFOIL en tiempo de ejecución ---
    results = [Xfoil.solve_alpha(a, re; iter=150, reinit=true) for a in alpha]
    local converged = [r[5] for r in results]
    
    c_l_clean = [converged[i] ? results[i][1] : NaN for i in 1:length(alpha)]
    c_d_clean = [converged[i] ? results[i][2] : NaN for i in 1:length(alpha)]
    c_m_clean = [converged[i] ? results[i][4] : NaN for i in 1:length(alpha)]
    
    for i in 2:length(alpha)
        if !isnan(c_l_clean[i]) && !isnan(c_l_clean[i-1]) && abs(c_l_clean[i] - c_l_clean[i-1]) > 0.15;  c_l_clean[i] = NaN; end
        if !isnan(c_d_clean[i]) && !isnan(c_d_clean[i-1]) && (c_d_clean[i] > 0.045 || abs(c_d_clean[i] - c_d_clean[i-1]) > 0.015); c_d_clean[i] = NaN; end
        if !isnan(c_m_clean[i]) && !isnan(c_m_clean[i-1]) && abs(c_m_clean[i] - c_m_clean[i-1]) > 0.015;  c_m_clean[i] = NaN; end
    end

    # Graficar curvas sólidas de tu simulación viva
    plot!(p1, alpha, c_l_clean, label="XFOIL ($lb)", color=clr, linewidth=1.5, legend_title="Reynolds Number")
    plot!(p2, alpha, c_d_clean, label="XFOIL ($lb)", color=clr, linewidth=1.5, legend_title="Reynolds Number")
    plot!(p3, alpha, c_m_clean, label="XFOIL ($lb)", color=clr, linewidth=1.5, legend_title="Reynolds Number")

    # --- PARTE B: Cargar, FILTRAR y graficar datos de Airfoil Tools (Líneas Punteadas) ---
    a_f, cl_f, cd_f, cm_f = leer_datos_archivo_limpios(nombre_txt)
    
    if !isempty(a_f)
        plot!(p1, a_f, cl_f, label="AirfoilTools ($lb)", color=clr, linestyle=:dash, linewidth=1.2)
        plot!(p2, a_f, cd_f, label="AirfoilTools ($lb)", color=clr, linestyle=:dash, linewidth=1.2)
        plot!(p3, a_f, cm_f, label="AirfoilTools ($lb)", color=clr, linestyle=:dash, linewidth=1.2)
    else
        println("Advertencia: No se encontró o está vacío el archivo $nombre_txt")
    end
end
=#
# 4. Renderizar de forma independiente cada ventana de gráfico
display(p1)
display(p2)
display(p3)

#=
#Use loops to get values for different thicknesses
#use a loop (for) to plot! the different thicknesses

#Generate the lift to drag ratio plot:
lift_over_drag = c_l ./ c_d 

p5 = plot(alpha, lift_over_drag,
    label="t = 10%",
    ylabel="Lift to Drag Ratio Cl/Cd",
    xlabel="Angle of Attack", ylims=(-110,115))
p6 = plot(alpha, c_l,
    label="t = 10%",
    ylabel="Lift coefficient Cl",
    xlabel="Angle of Attack", ylims=(-1.7,1.9))

#=create plots for lift to drag ratio and lift coefficient 
for camber chaning=#
p8 = plot(alpha, lift_over_drag,
    label="c = 1.0%",
    ylabel="Lift to Drag Ratio Cl/Cd",
    xlabel="angle of attack")

p9 = plot(alpha, c_l,
    label="c = 1.0%",
    ylabel="Lift coefficient Cl",
    xlabel="angle of attack")

# Suing 11% max thickness instead of 10%
function run_xfoil(x, y, re, c_t)
    xx = copy(x)
    yy = copy(y)
    for i in 1:div(length(yy), 2)
        y_bar      = (yy[i] + yy[end-i+1])/2
        t          = yy[i] - yy[end-i+1]
        yy[i]       = y_bar + (t*c_t)/2
        yy[end-i+1] = y_bar - (t*c_t)/2
    end 
    Xfoil.set_coordinates(xx, yy)
    xr, yr = Xfoil.pane()
    for i = 1:n_a
        c_l[i], c_d[i], c_dp[i], c_m[i], converged[i] = Xfoil.solve_alpha(alpha[i], 
        re; iter=100, reinit=true)  
    end
    return c_l, c_d, c_dp, c_m 
end


c_t = [(2.0, :pink, "t = 20%"), (3.0, :orange, "t = 30%"), (4.0, :red, "t = 40%")]

for (i, clr, lb) in c_t
    c_l, c_d, c_dp, c_m = run_xfoil(x, y, re, i)
    lift_over_drag = c_l ./ c_d
    plot!(p5, alpha, lift_over_drag,
        label=lb, color=clr)
    plot!(p6, alpha, c_l,
        label=lb, color=clr)
end

display(p5)
display(p6)

#now varying the camber instead of the thickness

function camber_scaling(x, y, re, c_c)
    xx = copy(x)
    yy = copy(y)
    #find the maximum y value for camber line
    camberline = zeros(length(yy))
    for i in 1:div(length(yy), 2)
        camberline[i] = (yy[i] + yy[end-i+1])/2
    end
    y_max = maximum(camberline)
    #Scale the camber line and adjust the upper and lower surfaces
    for i in 1:div(length(yy), 2)
        y_bar       = (yy[i] + yy[end-i+1])/2
        t           = yy[i] - yy[end-i+1]
        k           = c_c/y_max
        y_bar       = y_bar*k
        yy[i]       = y_bar + (t/2)
        yy[end-i+1] = y_bar - (t/2)
    end 
    Xfoil.set_coordinates(xx, yy)
    xr, yr = Xfoil.pane()
    for i = 1:n_a
        c_l[i], c_d[i], c_dp[i], c_m[i], converged[i] = Xfoil.solve_alpha(alpha[i], 
        re; iter=100, reinit=true)  
    end
    return c_l, c_d, c_dp, c_m 
end



#Create loop for plotting
c_c = [(0.003, :black, "c = 0.3%"), (0.005, :green, "c = 0.5%"), (0.04, :pink, "c = 4%"), (0.07, :orange, "c = 7%")]

for (camber_value, clr, lb) in c_c
    c_l, c_d, c_dp, c_m = camber_scaling(x, y, re, camber_value)
    lift_over_drag = c_l ./ c_d
    plot!(p8, alpha, lift_over_drag,
        label=lb, color=clr)
    plot!(p9, alpha, c_l,
        label=lb, color=clr)
end

display(p8)
display(p9) 
=#