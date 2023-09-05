using DataDrivenDiffEq
using OrdinaryDiffEq
using DataDrivenDMD
using LinearAlgebra
using Plots

# system to be used
mu = -0.05
λ = -1
A = [mu 0 0; 0  λ -λ; 0 0  2mu]

# initial conditions
u0A = [1.5; -1; 2.25]
u0B = [1; -1; 1]
u0C = [2; -1; 4]
tspan = (0.0, 1000.0)

sysA = ODEProblem((u,p,t) -> A*u, u0A, tspan)
sysB = ODEProblem((u,p,t) -> A*u, u0B, tspan)
sysC = ODEProblem((u,p,t) -> A*u, u0C, tspan)

solA = solve(sysA, Tsit5(), dt=0.1, saveat=0.1)
solB = solve(sysB, Tsit5(), dt=0.1, saveat=0.1)
solC = solve(sysC, Tsit5(), dt=0.1, saveat=0.1)

# start the data-driven model
X = Array(solA)
t = solA.t
ddprob = ContinuousDataDrivenProblem(X, t)

## look at koopman trajectories
# plot the solution
plot3d(solA[1,:], solA[2,:], solA[3,:], color = :red, label = "Solution A", camera = (-25, 30))
plot3d!(solB[1,:], solB[2,:], solB[3,:], color = :blue, label = "Solution B")
plot3d!(solC[1,:], solC[2,:], solC[3,:], color = :green, label = "Solution C")
