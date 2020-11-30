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

τ=1.0
function lotka_voltera!(du,u,h,p,t)
    🐰,🐺=u
    🕐🐰 = h(p,t-τ;idxs=1)
    α,β,γ,δ = p
    du[1] = d🐰 = α*🕐🐰-β*🐰*🐺
    du[2] = d🐺 = γ*🐰*🐺-δ*🐺
end

u₀ = [1.0,1.0]
tspan = (0.0,10.0)
h(p,t) = [1.0,1.0]
h(p,t;idxs=1)=1.0
p = [1.5,1.0,3.0,1.0]
probdel = DDEProblem(lotka_voltera!,u₀,h,tspan,p,constant_lag=[τ])
soldel = solve(probdel)
plot(soldel)

# using the callback functions
🔥🐰_condition(u,t,integrator) = u[2] - 4
🔥🐰_affect!(integrator) = integrator.u[2] -= 1
🔥🐰_cb = ContinuousCallback(🔥🐰_condition,🔥🐰_affect!)

solcall = solve(probdel,Tsit5(), callback = 🔥🐰_cb)
plot(solcall)


scatter(sol.t,dataset')
tmp_prob = remake(prob,p=[1.2,0.8,2.5,0.8])
# remake the solver
tmp_sol = solve(tmp_prob,Tsit5(),saveat=0.1)
plot(tmp_sol)
dataset = Array(tmp_sol)


function loss(p)
    tmp_prob = remake(prob,p=p)
    tmp_sol = solve(tmp_prob,saveat=0.1)
    sum(abs2,Array(tmp_sol)-dataset),tmp_sol
end

using DiffEqFlux, Optim
pinit = [1.2,0.8,2.5,0.8]
res = DiffEqFlux.sciml_train(loss,pinit,BFGS())

function plot_callback(p,l,tmp_sol)
    @show l
    tmp_prob = remake(prob,p=p)
    tmp_sol = solve(tmp_prob,saveat=0.1)
    fig = plot(tmp_sol)
    scatter!(fig,sol.t,dataset')
    display(fig)
    false
end
using Flux
pinit = [1.3,0.9,2.5,0.4]
res = DiffEqFlux.sciml_train(loss,pinit,ADAM(0.01),cb=plot_callback,maxiters=200)
res.minimizer

using BlackBoxOptim
res = DiffEqFlux.sciml_train(loss,pinit,DiffEqFlux.BBO(),lower_bounds=0.5ones(4),upper_bounds=4.0ones(4))
