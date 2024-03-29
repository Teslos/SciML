# pendulum implemented as Hamiltonian
using DiffEqPhysics
using DifferentialEquations
using Plots
default(fmt = :png, size = (450, 250))
using StaticArrays

module My
mutable struct ParamSinglePendulum{Tg, Tl}
    g::Tg
    l::Tl
end
end

function singlependulum(p, q, param, t = nothing)
    (; g, l) = param
    dθ = p[1]
    θ = q[1]
    (1/2)*(dθ/l)^2 - g*l*cos(θ)
end

plot_singlependulum(sol) = plot(sol; label=["dθ/dt" "θ"], legend=:outertopright)

param = My.ParamSinglePendulum(9.8, 1.0) # on the earth
q0 = SVector(0.0)
p0 = SVector(π/2)
tspan = (0.0, 10.0)
prob = HamiltonianProblem(singlependulum, p0, q0, tspan, param)

integrator = init(prob, Tsit5())
param.g = 9.8 # on the earth
dt = 0.1
for _ in 1:40
    step!(integrator, dt, true)
end
plot_singlependulum(integrator.sol)

param.g = 1.6 # on the moon
dt = 0.1
for _ in 1:210
    step!(integrator, dt, true)
end
plot_singlependulum(integrator.sol)

param.g = 9.8 # on the earth
dt = 0.1
for _ in 1:200
    step!(integrator, dt, true)
end
plot_singlependulum(integrator.sol)