# this file we try to solve Hamiltonian based on 2d hexagonal lattice
using DifferentialEquations
using DiffEqPhysics

using Plots

# initial conditions for the 2d hexagonal oscillators
p0 = [0.0 0.0 ; 
      0.0 0.0 ; 
      0.0 0.0 ; 
      0.0 0.0 ; 
      0.0 0.0 ; 
      0.0 0.0 ; 
      0.0 0.0 ]

q0 = [0.0 0.0 ; 
      0.0 1.0 ; 
      0.0 0.0 ; 
      0.0 0.0 ; 
      0.0 0.0 ; 
      0.0 0.0 ; 
      0.0 0.0 ]

# this tediously define the neigh. in the hexagonal mesh in 3-directions.
Σ = [0 0 0 0 0 0 1;
     0 0 0 0 0 0 1;
     0 0 0 0 0 0 1;
     0 0 0 0 0 0 1;
     0 0 0 0 0 0 1;
     0 0 0 0 0 0 1;
     1 1 1 1 1 1 0;]
# positions of the oscillators
positions = [1.0 0.0 ; 
             0.0 1.0 ; 
            -1.0 -1.0 ; 
            -1.0 0.0; 
             0.0 -1.0; 
             1.0 1.0 ; 
             0.0 0.0 ]

tspan = (0.0,100.0)

# Total energy of the system is consists of the potential energy
function V(σ₁, σ₂, σ₃)
    v = 0.0; α = 0.25
    for i = 1:length(σ₁)-1
        v += α/3 * (σ₁[i+1]-σ₁[i])^3
        v += 1/2 * (σ₁[i+1]-σ₁[i])^2
    end

    for j = 1:length(σ₂)-1
        v += α/3 * (σ₂[j+1]-σ₂[j])^3
        v += 1/2 * (σ₂[j+1]-σ₂[j])^2
    end

    for k = 1:length(σ₃)-1
        v += α/3 * (σ₃[k+1]-σ₃[k])^3
        v += 1/2 * (σ₃[k+1]-σ₃[k])^2
    end
    return v
end

function Vn(σ₁, σ₂)
    v = 0.0; α = 0.25
    n, m = size(Σ)
    σ₃ = -(σ₁ .+ σ₂)  # this is the third component of the vector
    for i=1:n
        for j=1:m
            if Σ[i,j] == 1
                v += α/3 * (σ₁[j]-σ₁[i])^3
                v += 1/2 * (σ₁[j]-σ₁[i])^2
                v += α/3 * (σ₂[j]-σ₂[i])^3
                v += 1/2 * (σ₂[j]-σ₂[i])^2
                v += α/3 * (σ₃[j]-σ₃[i])^3
                v += 1/2 * (σ₃[j]-σ₃[i])^2
            end
        end
    end
    return v
end

T(pσ₁, pσ₂) = 1//2 * sum(pσ₁ .^2 + pσ₂ .^2)
E(σ₁,σ₂,pσ₁,pσ₂) = T(pσ₁,pσ₂) + Vn(σ₁,σ₂)
# define the function for hexagonal lattice
function hexagonal2d(p,q,param,t)
    σ₁ = q[:,1]
    σ₂ = q[:,2]
    
    dσ₁ = p[:,1]
    dσ₂ = p[:,2]
    h = E(σ₁,σ₂,dσ₁,dσ₂)
    return h 
end

dt = 0.1
# Pass to solvers
prob = HamiltonianProblem(hexagonal2d, p0, q0, tspan, param = nothing)
integrator = init(prob, Tsit5(), abs_tol = 1e-16, rel_tol = 1e-16)

for _ in 1:1000
    step!(integrator, dt, true)
end


#plot the orbit
plot(integrator.sol.t, integrator.sol.u[1:2, :], title="The position of oscillators")

function hextoxy(h)
    origin = (0.0,0.0)
    y = h[1] - h[2] * sin(π/6) 
    x = - h[2] * cos(π/6)
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
        xy = hextoxy.(eachrow(sol.u[i].x[2]+positions))
        plot(xy, size=(400,300), xlim=(-axis_lim,axis_lim), ylim=(-axis_lim,axis_lim), markersize = 10, markershape = :circle,label ="",axis = [],
             title = str, title_location = :left)
     
        frame(anim)
    end
    gif(anim, fps = 30)
end

make_pretty_gif(integrator.sol)