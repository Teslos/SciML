using DifferentialEquations
using Plots
using LinearAlgebra
using FFTW
using DataDrivenDiffEq
using LinearAlgebra
using OrdinaryDiffEq
using DataDrivenDMD

# This simulation is trying to implement the logical switch
# using the driven duffing oscillator. Idea is to use the network to learn 
# the logical operation. Similar to implementation of the XOR gate using NN.

# these are the parameters used in original paper
# by Jothimurugan Thamilmaran et.al.
d = 0.1
ω0 = 1.0
β = 20
δ = 1
f = 0.1

# this is same as previous but the
# update is done using first order diff. equations.
# in this equation we have two oscillators which are coupled
# with linear springs.
function duffing_driven(du,u,p,t)
    # coupling is done with linear springs
    f1,ω1,f2,ω2,ϕ = p
    du[1] = u[2]
    du[2] = -d*u[2] - (ω0^2)*u[1] - δ*(u[1]-u[3]) -β*u[1]^3 + f1 * cos(ω1*t+ϕ)
    du[3] = u[4]
    du[4] = -d*u[4] - (ω0^2)*u[3] + δ*(u[1]-u[3]) -β*u[3]^3 + f2 * cos(ω2*t+ϕ)
end

u0 = [0.1, 0.0, 0, 0.0]
tspan = (0.0, 10.0)
Δt = 1e-2
t = tspan[1]:Δt:tspan[end]
ω = [1.28, 1.58]    # this is driving frequency
ωa = [1.43, 1.42]   # this is the anti-resonance frequency

# using adjoint method to find smallest driving frequency
# that gives maximum amplitude
p = [0,ωa[1],0,ωa[1],0]
prob = ODEProblem{true, SciMLBase.NoSpecialize}(duffing_driven, u0, tspan, p)
sol = solve(prob, Tsit5(), adaptive=false, dt=Δt)
plot(sol,vars=(1,2),color=:red)
plot!(sol,vars=(3,4),color=:blue)

ddprob = DataDrivenProblem(sol)
@parameters t
@variables u(t)[1:4]
Ψ = Basis([u; u[1]^2; u[1]^3; u[3]^2; u[3]^3], u, independent_variables = t)
res = solve(ddprob, Ψ, DMDSVD(),digits=1)

dof(res)
basis = get_basis(res)
p = get_parameter_map(basis)
odes = get_results(res)
# try to solve the problem with found basis
prob2 = ODEProblem(basis, u0, tspan)
# sol2 = solve(prob2, Tsit5(), dt=0.1, saveat=0.1)

t = Symbolics.unwrap(get_iv(basis))

# create the system of equations from the basis
eqs = map(equations(basis)) do eq  
    eq.lhs ~ eq.rhs
end

@named sys = ODESystem(eqs, get_iv(basis), states(basis), parameters(basis))

x0 = [u[1] => u0[1], u[2] => u0[2], u[3] => u0[3], u[4] => u0[4]]
ps = get_parameter_map(basis)

ode_prob = ODEProblem(sys, x0, tspan, ps)
estimate = solve(ode_prob, Tsit5(), dt=0.01, saveat= Δt)
scatter!(estimate, vars=(1,2), label="estimate", makerkersize = 1.0, color = :red)
scatter!(estimate, vars=(3,4), label="estimate", makerkersize = 1.0, color = :blue)