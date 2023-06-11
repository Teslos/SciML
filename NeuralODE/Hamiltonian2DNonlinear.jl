using DifferentialEquations
using DiffEqPhysics

using Plots

# initial conditions for the 2d lattice oscillators
p0 = [0.0 0.1 ; 0.0 0.0; 0. 0.0; 0. 0.0; 0. 0.0; 0. 0.0; 0. 0.0; 0. 0.0;0. 0.0]
q0=  [0.5 0; 0.0 0.0;  0  0;0 0; 0 0; 0 0;0 0;0 0;0 0]
positions = [0.0 0; 1.0 0.0; 2  0;0 1;1 1; 2 1;0 2;1 2;2 2]
tspan = (0.,100.)

# Total energy of the system is consists
# this is potential energy
function V(x, y)
    v = 0.0; alpha=0.25
    for i=1:length(x)-1
        v += alpha/3 * (x[i+1] - x[i])^3
        v += 1/2 * (x[i+1] - x[i])^2
    end
    for j=1:length(y)-1
        v += alpha/3 * (y[j+1] - y[j])^3
        v += 1/2 * (y[j+1] - y[j])^2
    end
    return v
end

T(px,py) = 1//2 * sum(px.^2 + py.^2)


E(x,y,px,py) = V(x,y)+T(px,py)

# define the function
function lattice2d(p,q,param,t)
    x = q[:,1]
    y = q[:,2]
    dx = p[:,1]
    dy = p[:,2]

    h = E(x,y,dx,dy)
    return h
end

dt = 0.1
#Pass to solvers
prob = HamiltonianProblem(lattice2d, p0, q0, tspan, param=nothing)
integrator = init(prob, Tsit5(), abs_tol = 1e-16, rel_tol = 1e-16)
for _ in 1:1000
    step!(integrator, dt, true)
end

# Plot orbit
plot(integrator.sol.t, integrator.sol.u[1:2,:], title="The position of oscillators", xaxis="x", yaxis="y", leg=false)
for i in 2:9 
    plot!(integrator.sol, vars=(i*2,i*2+1), leg=false)
end
plot(integrator.sol, vars=(6,7), title="Phase space for the first oscillator", xaxis="Position", yaxis="Velocity")

energy = map(x->E(x...), integrator.sol.u)
@show ΔE = energy[1]-energy[end]

anim = @animate for i = 1:10:1000
    plot(integrator.sol.u[i, :])
end

gif(anim, "oscillator2d.gif", fps=10)
plot(integrator.sol.t, energy .- energy[1], title="Change in the energy over time", xaxis= "Time in the iterations", yaxis="Change in the energy")

#Animation
#==========================================================#
function make_pretty_gif(sol)
    
    timepoints = sol.t
    x = Vector{Vector{Float64}}(undef, 9)
    y = Vector{Vector{Float64}}(undef, 9)
    print(typeof(sol[1,:]))
    for i=1:9
        x[i] = sol[2*i+1,:]
        y[i] = sol[2*i+2,:]
    end

    axis_lim = 3

    anim = Animation()
    for i =1:length(timepoints)
        str = string("Time = ", round(timepoints[i],digits=1), " sec")
        plot([x[1][i]+positions[1,1]], [y[1][i]+positions[1,2]], size=(400,300), xlim=(-axis_lim,axis_lim), ylim=(-axis_lim,axis_lim), markersize = 10, markershape = :circle,label ="",axis = [])
        for j = 2:9 
         plot!([x[j][i]+positions[j,1]], [y[j][i]+positions[j,2]], markersize = 10, markershape = :circle,label ="",title = str, title_location = :left)
        end
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

sol1 = integrator.sol
make_pretty_gif(sol1)