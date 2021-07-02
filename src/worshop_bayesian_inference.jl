using Turing, Distributions, DataFrames, DifferentialEquations, DiffEqSensitivity

using MCMCChains, Plots, StatsPlots

function lotka_voltera!(du,u,p,t)
    🐰,🐺=u
    α,β,γ,δ = p
    du[1] = d🐰 = α*🐰-β*🐰*🐺
    du[2] = d🐺 = γ*🐰*🐺-δ*🐺
end
u0 = [1.0,1.0]
p = [1.5,1.0,3.0,1.0]
tspan = (0.0,10.0)
prob = ODEProblem(lotka_voltera!,u0,tspan,p)
odedata = Array(solve(prob, Tsit5(),saveat=0.1))
Turing.setadbackend(:forwarddiff)
@model function fitlv(data)
    σ ~ InverseGamma(2,3)
    α ~ truncated(Normal(1.3,0.5),0.5,2.5)
    β ~ truncated(Normal(1.2,0.5),0,2)
    γ ~ truncated(Normal(2.7,0.5),0,2)
    δ ~ truncated(Normal(1.3,0.5),0,2)

    p = [α,β,γ,δ]
    prob = ODEProblem(lotka_voltera!,u0,(0.0,10.0),p)
    predicted = solve(prob,Tsit5(),saveat=0.1)

    for i = 1:length(predicted)
        data[:,i] ~ MvNormal(predicted[i],σ)
    end
end

model = fitlv(odedata)
chain = sample(model,NUTS(.36),10000)
plot(chain)
