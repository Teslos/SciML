using DiffEqFlux, Optimization, OptimizationFlux, DifferentialEquations, LinearAlgebra
k, α, β, γ = 1, 0.1, 0.2, 0.3

function dxdt_train(du,u,p,t)
    du[1] = u[2]
    du[2] = -k*u[1]-α*u[1]^3-β*u[2]-γ*u[2]^3
end

u0 = [1.0,0.0]
ts = collect(0.0:0.1:tspan[2])
prob_train = ODEProblem{true}(dxdt_train, u0, tspan)
data_train = Array(solve(prob_train, Tsit5(), saveat = ts))

A = [LegendreBasis(10), LegendreBasis(10)]
nn = TensorLayer(A,1)

f = x->min(30one(x),x)

function dxdt_pred(du,u,p,t)
    du[1] = u[2]
    du[2] = -p[1]*u[1]-p[2]*u[2] + f(nn(u,p[3:end])[1])
end

𝒜 = zeros(102)

prob_pred = ODEProblem{true}(dxdt_pred, u0, tspan)

function predict_adjoint(θ)
    x = Array(solve(prob_pred, Tsit5(), p = θ, saveat = ts,
    sensealg = InterpolatingAdjoint(autojacvec = ReverseDiffVJP(true))))
end

function loss_adjoint(θ)
    x = predict_adjoint(θ)
    loss = sum(norm.(x - data_train))
    return loss
end

iter = 0
function dcallback(θ, l)
    global iter
    iter += 1
    if iter % 10 == 0
        println(l)
    end
    return false
end

adtype = Optimization.AutoZygote()
optf = Optimization.OptimizationFunction((x,p) -> loss_adjoint(x), adtype)
optprob = Optimization.OptimizationProblem(optf, 𝒜)
res1 = Optimization.solve(optprob, ADAM(0.05), callback = dcallback, maxiters = 150)

optprob2 = Optimization.OptimizationProblem(optf, res1.u)
res2 = Optimization.solve(optprob2, ADAM(0.001), callback = dcallback, maxiters = 200)

using Plots
data_pred = predict_adjoint(res1.u)
plot(ts, data_train[1,:],label="X (ODE)")
plot!(ts, data_train[2,:],label="V (ODE)")
plot!(ts, data_pred[1,:], label="X (NN)")
plot!(ts, data_pred[2,:],label="V (NN)")