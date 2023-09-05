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

# define the nonlinear system
function nonlinear_system(du, u, p, t)
    du[1] = p[1]*u[1]
    du[2] = p[2]*u[2]- p[2]*u[1]^2
end

u0 = [1.5; -1]
tspan = (0.0, 1000.0)
p = [-0.05; -1.0]

problem = ODEProblem{true, SciMLBase.NoSpecialize}(nonlinear_system, u0, tspan, p)
sol = solve(problem, Tsit5(), dt=0.1, saveat=0.1)
plot!(sol)

prob = DataDrivenProblem(sol)
@parameters t
@variables u(t)[1:2]
Ψ = Basis([u; u[1]^2], u, independent_variables = t)
res = solve(prob, Ψ, DMDPINV(),digits=2)

dof(res)
basis = get_basis(res)
get_parameter_map(basis)

plot(res)