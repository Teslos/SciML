# Extended Dynamic Mode Decomposition
using DataDrivenDiffEq
using LinearAlgebra
using OrdinaryDiffEq
using DataDrivenDMD
using Plots

# integrate the Koopman trajectories
function slow_manifold(du, u, p, t)
    du[1] = p[1]*u[1]
    du[2] = p[2]*u[2] - p[2]*u[3]
    du[3] = 2*p[1]*u[3]
end

u0 = [1.5; -1; 2.25]
u1 = [1; -1; 1]
u2 = [2; -1; 4]

tspan = (0.0, 1000.0)

p = [-0.05; -1.0]   

problem = ODEProblem{true, SciMLBase.NoSpecialize}(slow_manifold, u0, tspan, p)
sol = solve(problem, Tsit5(), dt=0.1, saveat=0.1)
plot(sol)

prob = DataDrivenProblem(sol)
plot(prob)

# the parameters are defined in the paper
δ=0.5
β=-1
α=1
p = [δ; β; α] 
# define the nonlinear Duffing system
function duffing(du, u, p, t)
    du[1] = u[2]
    du[2] = -p[1]*u[2]- p[2] * u[1] - p[3] * u[1]^3 
end

# define the initial conditions for x coordinate and v velocity between -2 and 2
using Distributions
x0 = Uniform(-2, 2)
v0 = Uniform(-2, 2)
u0 = [rand(x0); rand(v0)]
tspan = (0.0, 10.0)

problem = ODEProblem{true, SciMLBase.NoSpecialize}(duffing, u0, tspan, p)
sol = solve(problem, Tsit5(), dt=0.1, saveat=0.1)
plot(sol,vars=(1,2))

prob = DataDrivenProblem(sol)
@parameters t
@variables u(t)[1:2]
Ψ = Basis([u; u[1]^2; u[1]^3], u, independent_variables = t)
res = solve(prob, Ψ, DMDSVD(),digits=1)

dof(res)
basis = get_basis(res)
p = get_parameter_map(basis)
odes = get_results(res)
# try to solve the problem with found basis
prob2 = ODEProblem(basis, u0, tspan, p)
sol2 = solve(prob2, Tsit5(), dt=0.1, saveat=0.1)