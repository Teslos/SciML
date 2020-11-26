using DifferentialEquations
using Plots
function lotka_voltera(du,u,p,t)
    🐰,🐺=u
    α,β,γ,δ = p
    du[1] = d🐰 = α*🐰-β*🐰*🐺
    du[2] = d🐺 = γ*🐰*🐺-δ*🐺
end
u₀=[1.0,1.0]
tspan = [0.0,10.0]
p = [1.5,1.0,3.0,1.0]
prob = ODEProblem(lotka_voltera,u₀,tspan, p)
sol = solve(prob,Tsit5(),abstol=1e-8,reltol=1e-8)

plot(sol)
plot(sol,vars=(1,2))

function multiplicative_noise!(du,u,p,t)
    🐰, 🐺 = u
    du[1] = 0.3*🐰
    du[2] = 0.3*🐺
end

prob = SDEProblem(lotka_voltera,multiplicative_noise!,u₀,tspan, p)
solp = solve(prob)
plot(solp)

ensembleprob = EnsembleProblem(prob)
solens = solve(ensembleprob,SOSRI(),EnsembleThreads(),trajectories=1000)
plot(solens)

summ = EnsembleSummary(solens)
plot(summ)
