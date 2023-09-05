using DataDrivenDiffEq
using OrdinaryDiffEq
using DataDrivenDMD
using Plots
using LinearAlgebra

A = [-0.9 0.2; 0.0 -0.2]
u0 = [10.0; -10.0]
tspan = (0.0, 10.0)

f(u, p, t) = A*u
sys = ODEProblem(f, u0, tspan)
sol = solve(sys, Tsit5(), saveat=0.05)

X = Array(sol)
t = sol.t
prob = ContinuousDataDrivenProblem(X, t)

plot(prob)

# look at the derivative of the solution
DX = Array(sol(t, Val{1}))
scatter(t, DX',  label=["Solution" nothing], color = :red, legend = :bottomright)
plot!(t, prob.DX', label=["Linear interpolation" nothing], color = :blue)

# do the DMD analysis
res = solve(prob, DMDSVD(), digits=1)   
println(get_basis(res))
plot(res)