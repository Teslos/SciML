# hexagonal lattice
using OrdinaryDiffEq, Plots
using ComponentArrays

# Parameters
α = 1.0
β = 0.0
γ = 1.0
d = 0.0
ω = 1.0
tspan = (0.0, 8π)
t = range(tspan[1], tspan[2], length = 200)
# initial conditions
u0 =[0.2 0.0 0.0; 0.2 0.2 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0  0.0; 0.0 0.0 0.0]
du0 = [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0  0.0; 0.0 0.0 0.0]
# positions of the masses (σ1, σ2, σ3)
x0 = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0; -1.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 -1.0; 0.0 0.0 0.0]

uc0 =[0.0 0.1; 0.0 0.0 ; 0.0 0.0 ; 0.0 0.0 ; 0.0 0.0 ; 0.0 0.0 ; 0.0 0.0]
duc0 = [0.0 0.0 ; 0.0 0.0 ; 0.0 0.0 ; 0.0 0.0 ; 0.0 0.0 ; 0.0 0.0  ; 0.0 0.0 ]

# positions of the masses (x,y)
xc0 = [0.0 1.0; -√3/2 -0.5 ; √3/2 -0.5 ; 0 -1.0; √3/2 0.5 ;  -√3/2 0.5 ; 0.0 0.0 ]
#define the problem
function hexagonal_duffing(ddu, du, u, param, t)
    α, β, d, ω = param
    ddu[1,1] = -α * (u[1,1] - u[7,1]) -β*(u[1,1] - u[7,1])^3 + d*cos(ω * t)
    ddu[2,2] = -α * (u[2,2] - u[7,2]) -β*(u[2,2] - u[7,2])^3 + d*cos(ω * t)
    ddu[3,3] = -α * (u[3,3] - u[7,3]) -β*(u[3,3] - u[7,3])^3 + d*cos(ω * t)
    ddu[4,1] = -α * (u[4,1] - u[7,1]) -β*(u[4,1] - u[7,1])^3 + d*cos(ω * t)
    ddu[5,2] = -α * (u[5,2] - u[7,2]) -β*(u[5,2] - u[7,2])^3 + d*cos(ω * t)
    ddu[6,3] = -α * (u[6,3] - u[7,3]) -β*(u[6] - u[7])^3 + d*cos(ω * t)

    ddu[7,1] = -α * (u[1,1] - u[7,1]) -α*(u[7,1] - u[4,1])
    ddu[7,2] = -α * (u[2,2] - u[7,2]) -α*(u[7,2] - u[5,2])
    ddu[7,3] = -α * (u[3,3] - u[7,3]) -α*(u[7,3] - u[6,3])

    ddu[7,1] += -β * (u[1,1] - u[7,1])^3 -β*(u[7,1] - u[4,1])^3 + d*cos(ω * t)
    ddu[7,2] += -β * (u[2,2] - u[7,2])^3 -β*(u[7,2] - u[5,2])^3 + d*cos(ω * t)
    ddu[7,3] += -β * (u[3,3] - u[7,3])^3 -β*(u[7,3] - u[6,3])^3 + d*cos(ω * t)
             
    return ddu
end

# define hexagonal lattice in cartesian coordinates
# all equations are explicitly written out.
function hexagonal_duffing_cart(ddu, du, u, param, t)
    α, β, d, ω = param
   
    #ddu[1,1] = -α * (u[1,1] - u[7,1]) -β*(u[1,1]-u[7,1])^3 + d*cos(ω * t)
    ddu[1,1] = 0.0
    ddu[1,2] = -α * (u[1,2] - u[7,2]) -β*(u[1,2]-u[7,2])^3 + d*cos(ω * t)

    ddu[2,1] = √3/2*(-α * (u[2,1] - u[7,2]) -β*(u[2,1]-u[7,1])^3) + d*cos(ω * t)
    ddu[2,2] =  1/2*(-α * (u[2,2] - u[7,2]) -β*(u[2,2]-u[7,2])^3) + d*cos(ω * t)

    ddu[3,1] = -√3/2*(-α * (u[3,1] - u[7,1]) -β*(u[3,1]-u[7,1])^3) + d*cos(ω * t)
    ddu[3,2] =   1/2*(-α * (u[3,2] - u[7,2]) -β*(u[3,2]-u[7,2])^3) + d*cos(ω * t)

    #ddu[4,1] = -α * (u[4,1] - u[7,1]) -β*(u[4,1]-u[7,1])^3 + d*cos(ω * t)
    ddu[4,1] = 0.0
    ddu[4,2] = -α * (u[4,2] - u[7,2]) -β*(u[4,2]-u[7,2])^3 + d*cos(ω * t)

    ddu[5,1] = -√3/2*(-α * (u[5,1] - u[7,1]) -β*(u[5,1]-u[7,1])^3) + d*cos(ω * t)
    ddu[5,2] =  -1/2*(-α * (u[5,2] - u[7,2]) -β*(u[5,2]-u[7,2])^3) + d*cos(ω * t)

    ddu[6,1] = √3/2*(-α * (u[6,1] - u[7,1]) -β*(u[6,1]-u[7,1])^3) + d*cos(ω * t)
    ddu[6,2] = -1/2*(-α * (u[6,2] - u[7,2]) -β*(u[6,2]-u[7,2])^3) + d*cos(ω * t)

    ddu[7,1] = #= -α * (u[7,1] - u[1,1]) -α * (u[7,1] - u[4,1]) =#  -α * -√3/2 * (u[7,1] - u[2,1])
               -α * √3/2*(u[7,1] - u[5,1]) -α * √3/2* (u[7,1] - u[3,1]) -α * -√3/2*(u[7,1] - u[6,1])
    ddu[7,1] += #= -β * (u[7,1] - u[1,1])^3 -β * (u[7,1] - u[4,1])^3 =# -β * -√3/2 * (u[7,1] - u[2,1])^3
                -β * √3/2*(u[7,1] - u[5,1])^3 -β * √3/2 * (u[7,1] - u[3,1])^3 -β * -√3/2* (u[7,1] - u[6,1])^3   
                +d*cos(ω * t)        
    ddu[7,2] = -α * (u[7,2] - u[1,2]) -α * (u[7,2] - u[4,2]) -α * -1/2*(u[7,2] - u[2,2])
               -α * 1/2*(u[7,2] - u[5,2]) -α * -1/2* (u[7,2] - u[3,2]) -α * 1/2*(u[7,2] - u[6,2])
    ddu[7,2] += -β * (u[7,2] - u[1,2])^3 -β * (u[7,2] - u[4,2])^3 -β * -1/2*(u[7,2] - u[2,2])^3
                -β * 1/2*(u[7,2] - u[5,2])^3 -β * -1/2*(u[7,2] - u[3,2])^3 -β * 1/2*(u[7,2] - u[6,2])^3   
                 +d*cos(ω * t)    

end

#pass to solvers
p = (α, β, d, ω)
prob = SecondOrderODEProblem(hexagonal_duffing_cart, duc0, uc0, tspan, p)
sol = solve(prob, Tsit5(), saveat = t)

struct Orientation
    f0::Float64
    f1::Float64
    f2::Float64
    f3::Float64
    b0::Float64
    b1::Float64
    b2::Float64
    b3::Float64
    startAngle::Float64
end

# assume pointy layout of the hexagon
function hexToPixel(h)
    origin = (0.0, 0.0)
    f0 = sqrt(3.0); f1 = sqrt(3.0) / 2.0
    f2 = 0.0; f3 = 3/2
    size = (sqrt(3.0)/2.0, 1.0)

    x = (f0 * h[1] + f1 * h[2]) * size[1]
    y = (f2 * h[1] + f3 * h[2]) * size[2]
    return x + origin[1], y + origin[2]
end

function hextoxy(h)
    origin = (0.0,0.0)
    y = h[1] - h[2] * cos(π/3) - h[3] * cos(π/3)
    x = h[3] * cos(π/6) - h[2] * cos(π/6)
    return (x,y)
end
function make_pretty_gif(sol)
    
    timepoints = sol.t
    v = Vector{Matrix{Float64}}(undef, 20)

    for i = 1:20
        v[i] = sol.u[i].x[2]
    end
    print("size:",size(v))
    x, y = hextoxy.(v)
    println("xy coordinates:")
    print(x)
    axis_lim = 3

    anim = Animation()
    for i =1:length(timepoints)
        str = string("Time = ", round(timepoints[i],digits=1), " sec")
        # convert triangular system to x, y coordinates
        #xy = hextoxy.(eachrow(sol.u[i].x[2]+x0))
        xy = sol.u[i].x[2]+xc0
        println(maximum(xy))
        plot(xy[:,1],xy[:,2], size=(400,300), xlim=(-axis_lim,axis_lim), ylim=(-axis_lim,axis_lim), markersize = 10, 
            markershape = :circle,label ="",axis = [], title = str, title_location = :left)
        # for j = 2:7
        #  x,y = hextoxy(v)    
        #  plot!([x], [y], markersize = 10, markershape = :circle,label ="",title = str, title_location = :left)
        # end
        #= if i > 8 #rainbow trail
            plot!([x2[i-2:i]],   [y2[i-2:i]],  alpha = 0.15, linewidth = 2, color = :red, label=nothing)
            plot!([x2[i-3:i-2]], [y2[i-3:i-2]],alpha = 0.15, linewidth = 2, color = :orange, label=nothing)
            plot!([x2[i-4:i-3]], [y2[i-4:i-3]],alpha = 0.15, linewidth = 2, color = :yellow, label=nothing)
            plot!([x2[i-6:i-4]], [y2[i-6:i-4]],alpha = 0.15, linewidth = 2, color = :green, label=nothing)
            plot!([x2[i-7:i-6]], [y2[i-7:i-6]],alpha = 0.15, linewidth = 2, color = :blue, label=nothing)
            plot!([x2[i-8:i-7]], [y2[i-8:i-7]],alpha = 0.15, linewidth = 2, color = :purple, label=nothing)
        end =#
        frame(anim)
    end
    gif(anim, fps = 30)
end

make_pretty_gif(sol)

