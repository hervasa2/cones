

magnitude(track) = sqrt(Transpose(track)*track)
export magnitude

function plot_point(point; markersize = 2)
    scatter([point[1]],[point[2]],[point[3]], markersize = markersize)
end
export plot_point

function plot_point!(point; markersize = 2)
    scatter!([point[1]],[point[2]],[point[3]], markersize = markersize)
end
export plot_point!

function generate_cone(axis, vertex, opening_angle, height; 
        draw = false, style = "parametric", palette = :magma, poly = 32, 
        height_mode = "relative", rand_cone = false, 
        circles = 500, start = 0, phase = 0, colors = 25, phaseZaxis = 0, phaseYaxis = 0, vertexsize = 2)
    
    if height_mode == "relative"
        height = height*magnitude(axis)
    end
    if rand_cone
        axis = [rand(Uniform(-1,1)), rand(Uniform(-1,1)),rand(Uniform(-1,1))]
        vertex = [0,0,0]
        height = 1
        poly = rand(append!(collect(2:6), collect(64:64)))
        circles = 300 
        phase = rand(Normal(0,pi/100)) 
        phaseYaxis = rand(Normal(0,pi/100)) 
        phaseZaxis = rand(Normal(0,pi/100))  
        palette = rand([:viridis,:magma,:ice,:phase])
    end
    step = (height-start)/circles
    lines = height/step + 1
    
    z_basis = [0,0,1]
    x_basis = [1,0,0]
    Rz(ϕ) = [cos(ϕ) -sin(ϕ) 0; sin(ϕ) cos(ϕ) 0; 0 0 1]
    Ry(ϕ) = [cos(ϕ) 0 sin(ϕ); 0 1 0; -sin(ϕ) 0 cos(ϕ)]
    Project_xy = [1 0 0; 0 1 0; 0 0 0]
    
    if draw plot(legend=false) end
    icirc = 0;
    r_max = abs(tan(opening_angle)*height)
    x_cone = zeros(0)
    y_cone = zeros(0)
    z_cone = zeros(0)
    for circ in start:step:height
        
        axis = Rz(phaseZaxis)*Ry(phaseYaxis)*axis
        θ_z = acos(Transpose(z_basis)*axis/(1*magnitude(axis)))
        ϕ = acos(Transpose(x_basis)*(Project_xy*axis)/(1*magnitude(Project_xy*axis)))
        if θ_z == 0 || (θ_z >= pi-0.00001 && θ_z <= pi+0.00001)
            ϕ = 0
        end
        if axis[2] > 0
            R = Rz(ϕ)*Ry(θ_z)
        elseif axis[2] < 0
            R = Rz(2*pi-ϕ)*Ry(θ_z) 
        else
            R = Rz(ϕ)*Ry(θ_z)              
        end
        
        if opening_angle >= pi/2
            circle_center = vertex + circ*axis/magnitude(axis)
            α = pi-opening_angle
        else
            circle_center = vertex - circ*axis/magnitude(axis)
            α = opening_angle
        end
        r = magnitude(circle_center .- vertex)*tan(α)
        
        x = zeros(0)
        y = zeros(0)
        z = zeros(0)
        if(style == "cartesian")
            xstep = 2*r_max/100
            append!(x,collect(-r:xstep:r))
            append!(y,sqrt.(r^2 .- x.^2))
            append!(z,zeros(2*length(x)))
            append!(x,x)
            append!(y,-y)
        end
        accumulated_phase = phase*icirc
        if(style == "parametric")
            β = collect(0+accumulated_phase:pi/(poly/2):2*pi+accumulated_phase)
            append!(x,r*cos.(β))
            append!(y,r*sin.(β))
            append!(z,zeros(length(x)))
        end
        
        x_t = similar(x)
        y_t = similar(y)
        z_t = similar(z)
        for i in 1:1:length(x)
            point = [x[i],y[i],z[i]]
            transformed_point = R*point+circle_center
            x_t[i] = transformed_point[1]
            y_t[i] = transformed_point[2]
            z_t[i] = transformed_point[3]
            append!(x_cone, x_t[i])
            append!(y_cone, y_t[i])
            append!(z_cone, z_t[i])
        end
        if(style == "parametric" && draw )
            plot!(x_t ,y_t ,z_t, color = cgrad(palette)[convert(Int64, floor(icirc/(lines/colors)))+1])
        end
        if(style == "cartesian" && draw )
            scatter!(x_t ,y_t ,z_t, markersize = 0.1, color = "black")
        end
        icirc = icirc + 1;
        #check orthogonality
        #println(axis[1]*(x_t[1]-circle_center[1]) + axis[2]*(y_t[1]-circle_center[2]) 
                    #+ axis[3]*(z_t[1]-circle_center[3]))
    end
    if !draw return x_cone, y_cone, z_cone end
    if draw plot_point!(vertex, markersize = vertexsize) end
end
export generate_cone            

function generate_cone!(axis, vertex, opening_angle, height; 
        draw = "true", style = "parametric", palette = :magma, poly = 32, 
        height_mode = "relative", rand_cone = false, 
        circles = 500, start = 0, phase = 0, colors = 20, phaseZaxis = 0, phaseYaxis = 0, vertexsize = 2)
                
    if height_mode == "relative"
        height = height*magnitude(axis)
    end
    if rand_cone
        axis = [rand(Uniform(-1,1)), rand(Uniform(-1,1)),rand(Uniform(-1,1))]
        vertex = [0,0,0]
        height = 1
        poly = rand(append!(collect(2:6), collect(64:64)))
        circles = 300 
        phase = rand(Normal(0,pi/100)) 
        phaseYaxis = rand(Normal(0,pi/100)) 
        phaseZaxis = rand(Normal(0,pi/100))  
        palette = rand([:viridis,:magma,:ice,:phase])
    end
    step = (height-start)/circles
    lines = height/step + 1
    
    z_basis = [0,0,1]
    x_basis = [1,0,0]
    Rz(ϕ) = [cos(ϕ) -sin(ϕ) 0; sin(ϕ) cos(ϕ) 0; 0 0 1]
    Ry(ϕ) = [cos(ϕ) 0 sin(ϕ); 0 1 0; -sin(ϕ) 0 cos(ϕ)]
    Project_xy = [1 0 0; 0 1 0; 0 0 0]

    icirc = 0;
    r_max = abs(tan(opening_angle)*height)
    for circ in start:step:height
        
        axis = Rz(phaseZaxis)*Ry(phaseYaxis)*axis
        θ_z = acos(Transpose(z_basis)*axis/(1*magnitude(axis)))
        ϕ = acos(Transpose(x_basis)*(Project_xy*axis)/(1*magnitude(Project_xy*axis)))
        if θ_z == 0 || (θ_z >= pi-0.00001 && θ_z <= pi+0.00001)
            ϕ = 0
        end
        if axis[2] > 0
            R = Rz(ϕ)*Ry(θ_z)
        elseif axis[2] < 0
            R = Rz(2*pi-ϕ)*Ry(θ_z) 
        else
            R = Rz(ϕ)*Ry(θ_z)              
        end
        
        if opening_angle >= pi/2
            circle_center = vertex + circ*axis/magnitude(axis)
            α = pi-opening_angle
        else
            circle_center = vertex - circ*axis/magnitude(axis)
            α = opening_angle
        end
        r = magnitude(circle_center .- vertex)*tan(α)
        
        x = zeros(0)
        y = zeros(0)
        z = zeros(0)
        if(style == "cartesian")
            xstep = 2*r_max/100
            append!(x,collect(-r:xstep:r))
            append!(y,sqrt.(r^2 .- x.^2))
            append!(z,zeros(2*length(x)))
            append!(x,x)
            append!(y,-y)
        end
        accumulated_phase = phase*icirc
        if(style == "parametric")
            β = collect(0+accumulated_phase:pi/(poly/2):2*pi+accumulated_phase)
            append!(x,r*cos.(β))
            append!(y,r*sin.(β))
            append!(z,zeros(length(x)))
        end
        
        x_t = similar(x)
        y_t = similar(y)
        z_t = similar(z)
        for i in 1:1:length(x)
            point = [x[i],y[i],z[i]]
            transformed_point = R*point+circle_center
            x_t[i] = transformed_point[1]
            y_t[i] = transformed_point[2]
            z_t[i] = transformed_point[3]
        end
        if(style == "parametric" && draw )
            plot!(x_t ,y_t ,z_t, color = cgrad(palette)[convert(Int64, floor(icirc/(lines/colors)))+1])
        end
        if(style == "cartesian" && draw )
            scatter!(x_t ,y_t ,z_t, markersize = 0.1, color = "black")
        end
        icirc = icirc + 1;
        #check orthogonality
        #println(axis[1]*(x_t[1]-circle_center[1]) + axis[2]*(y_t[1]-circle_center[2]) 
                    #+ axis[3]*(z_t[1]-circle_center[3]))
    end
    if draw plot_point!(vertex, markersize = vertexsize) end
end
export generate_cone!