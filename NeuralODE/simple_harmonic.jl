# Simple Harmonic Oscillator problem
using OrdinaryDiffEq, Plots
using ComponentArrays, DiffEqFlux, Lux, Zygote, Random, Optimization, RecursiveArrayTools

#Parameters
ω = 1

# Initial Conditions
x0 = [0.0]
dx0 = [π/2]
tspan = (0.0, 2π)
t = range(tspan[1], tspan[2], length = 20)

ϕ = atan((dx0[1]/ω)/x0[1])

A = sqrt(x0[1]^2 + dx0[1]/2)

# Define the problem
function harmonic(ddu,du, u, ω, t)
    ddu .= -ω^2 * u
end

# Pass to solvers
prob = SecondOrderODEProblem(harmonic, dx0, x0, tspan, ω)
sol = solve(prob, DPRKN6(),saveat=t)
ode_data = Array(sol)
plot(sol, vars=[2,1], linewidth=2, title="Simple Harmonic Oscillator", xaxis = "Time", yaxis = "Elongation", label = ["x", "dx"])

rng = Random.default_rng()

u0 = [0.0]; du0=[ π/2]
model = Lux.Chain(Lux.Dense(1,50,tanh), Lux.Dense(50,1))
p, st = Lux.setup(rng, model)
p =  ComponentArray(p)
ff(du,u,p,t) = model(u,p,st)[1]
prob2 = SecondOrderODEProblem{false}(ff,du0, u0,tspan, p)

function prediction(p)
    Array(solve(prob2, Tsit5(), p=p, saveat = t, sensealg = InterpolatingAdjoint(autojacvec = ZygoteVJP())))
end

function loss_n_pde(p)
    pred = prediction(p)
    loss = sum(abs2, ode_data .- pred )
    return loss, pred
end
callback = function(p,l,pred)
    @show l
    l < 0.01 && Lux.stop()
end

using OptimizationOptimisers
optfunc = Optimization.OptimizationFunction((x, p) -> loss_n_pde(x), Optimization.AutoZygote() )
optprob = Optimization.OptimizationProblem(optfunc, p)
res = Optimization.solve(optprob, OptimizationOptimisers.Adam(0.01), callback = callback, maxiters = 100 )
l2 = loss_n_pde(res.minimizer)
