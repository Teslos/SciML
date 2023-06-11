using DifferentialEquations
using DiffEqPhysics

using Plots

p0 = [0.,0.1]; q0=[0.5,0]
tspan = (0.,100.)

# Total energy of the system is consists
V(x,y) = 1//2 * (x^2 + y^2 + 2x^2*y - 2//3 * y^3)
E(x,y,dx,dy) = V(x,y)+1//2*(dx^2 + dy^2)

# define the function
function Henon_Heiles(p,q,param,t)
    x = q[1]
    y = q[2]
    dx = p[1]
    dy = p[2]

    h = E(x,y,dx,dy)
    return h
end

dt = 0.1
#Pass to solvers
prob = HamiltonianProblem(Henon_Heiles, p0, q0, tspan, param)
integrator = init(prob, Tsit5(), abs_tol = 1e-16, rel_tol = 1e-16)
for _ in 1:1000
    step!(integrator, dt, true)
end

# Plot orbit
plot(integrator.sol, vars=(1,2), title="The orbit of the Henon Heiles system", xaxis="x", yaxis="y", leg=false)
plot(integrator.sol, vars=(1,3), title="Phase space for the Henon Heiles system", xaxis="Position", yaxis="Velocity")

energy = map(x->E(x...), integrator.sol.u)
@show ΔE = energy[1]-energy[end]

plot(integrator.sol.t, energy .- energy[1], title="Change in the energy over time", xaxis= "Time in the iterations", yaxis="Change in the energy")