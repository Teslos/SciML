# Simple Harmonic Oscillator problem
using OrdinaryDiffEq, Plots
using ComponentArrays, DiffEqFlux, Lux, Zygote, Random, Optimization, RecursiveArrayTools

#Parameters
const g = 9.81
L = 1.0

# Initial Conditions

u0 = [0, π/2]
tspan = (0.0, 2π)
t = range(tspan[1], tspan[2], length = 20)

ϕ = atan((dx0[1]/ω)/x0[1])

A = sqrt(x0[1]^2 + dx0[1]/2)

# Define the problem
function harmonic_nonlinear(du, u, p, t)
    θ = u[1]
    dθ = u[2]
    du[1] = dθ
    du[2] = -(g/L) * sin(θ)
end

# Pass to solvers
prob = ODEProblem(harmonic_nonlinear, u0, tspan, p)
sol = solve(prob, Tsit5(),saveat=t)
ode_data = Array(sol)
plot(sol, vars=[2,1], linewidth=2, title="Nonlinear Harmonic Oscillator", xaxis = "Time", yaxis = "Elongation", label = ["x", "dx"])

rng = Random.default_rng()

model = Lux.Chain(Lux.Dense(2,50,tanh), Lux.Dense(50,2))
p, st = Lux.setup(rng, model)
p =  ComponentArray(p)
ff(u,p,t) = model(u,p,st)[1]
prob2 = ODEProblem{false}(ff, u0,tspan, p)

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
res = Optimization.solve(optprob, OptimizationOptimisers.Adam(0.01), callback = callback, maxiters = 300 )


l2 = loss_n_pde(res.minimizer)
plot(l2[2]')
