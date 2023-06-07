using DiffEqPhysics
using DifferentialEquations
using Plots
default(fmt = :png, size = (450, 250))
using StaticArrays

# constants and setup
const m₁, m₂, L₁, L₂ = 1,2,1,2
const g = 9.81
initial = [0, π/3, 0, 3π/5]
tspan = (0.,50.)

# Convenvience function for transforming from polar to Cartesian coordinates
function polar2cart(sol; dt=0.02, l1 = L₁, l2 = L₂, vars=(2,4))
    u =  sol.t[1]:dt:sol.t[end]

    p1 = l1*map(x->x[vars[1]], sol.(u))
    p2 = l2*map(x->x[vars[2]], sol.(u))

    x1 = l1*sin.(p1)
    y1 = l1*-cos.(p1)

    (u, (x1 + l2*sin.(p2),y1-l2*cos.(p2)))
end

# define the double pendulum Hamiltonian.
# https://diego.assencio.com/?index=e5ac36fcb129ce95a61f8e8ce0572dbf
function double_pendulum_hamiltonian(p, q, param, t = nothing)
    θ = q
    h = (m₂*L₂^2*p[1]^2 + (m₁ + m₂)*L₁^2*p[2]^2 - 2*m₂*L₁*L₂*p[1]*p[2]*cos(θ[1]-θ[2]))/
    (2*m₂*L₁^2*L₂^2*(m₁+m₂*sin(θ[1]-θ[2])^2)) - (m₁ + m₂)*g*L₁*cos(θ[1]) - m₂*g*L₂*cos(θ[2])
    return h
end


p0 = [0.0,0.0]
q0 = [π/3,3π/5]
prob = HamiltonianProblem(double_pendulum_hamiltonian, p0, q0, tspan, param)

integrator = init(prob, Tsit5())
dt = 0.01
for _ in 1:4000
    step!(integrator, dt, true)
end
ts,ps = polar2cart(integrator.sol, l1=L₁,l2=L₂,dt=0.01)
plot(ps...)