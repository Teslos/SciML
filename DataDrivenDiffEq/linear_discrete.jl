using DataDrivenDiffEq
using OrdinaryDiffEq
using DataDrivenDMD
using Plots
using LinearAlgebra

A = [0.9 -0.2; 0.0 0.2]
u0 = [10.0; -10.0]
tspan = (0.0, 10.0)

f(u,p,t) = A*u

sys = DiscreteProblem(f, u0, tspan)
sol = solve(sys, FunctionMap())

prob = DataDrivenProblem(sol)

res = solve(prob, DMDSVD(), digits=1)
get_basis(res)

println(get_basis(ddsol))

plot(res)
plot(sol, label = string.([:x1, :x2]))
scatter!(prob)