using OrdinaryDiffEq, Plots
f1 = function(du,u,p,t)
    du[1] = p[1]*u[1]-p[2]*u[1]*u[2]
    du[2] = -p[3]*u[2]+p[4]*u[1]*u[2]
end
θ = [1.5,1.0,3.0,1.0]
u0 = [1.0;1.0]
tspan = (0.0,10.0)
prob1 = ODEProblem(f1,u0,tspan,θ)
sol = solve(prob1,Tsit5())
plot(sol)

using Distributions

θ = [Uniform(0.5,1.5),Beta(5,1),Normal(3,0.5),Gamma(5,2)]

_θ = rand.(θ)
prob1 = ODEProblem(f1,u0,tspan,_θ)
sol = solve(prob1,Tsit5())
plot(sol)

prob_func = function(prob,i,repeat)
    remake(prob,p=rand.(θ))
end

ensemble_prob = EnsembleProblem(ODEProblem(f1,u0,tspan,θ),
    prob_func=prob_func)
sol = solve(ensemble_prob, Tsit5(), EnsembleThreads(),trajectories=1000)
using DiffEqBase.EnsembleAnalysis
plot(EnsembleSummary(sol))
