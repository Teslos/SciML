using DiffEqFlux, OrdinaryDiffEq, Flux, Optim, Plots, DiffEqSensitivity

u0 = Float32[2.0; 0.0]
datasize = 30
tspan = (0.0f0, 1.5f0)
tsteps = range(tspan[1], tspan[2], length = datasize)

function trueODEfunc(du, u, p, t)
    true_A = [-0.1 2.0; -2.0 -0.1]
    du .= ((u.^3)'true_A)'
end

prob_trueode = ODEProblem(trueODEfunc, u0, tspan)
ode_data = Array(solve(prob_trueode, Tsit5(), saveat = tsteps))

dudt2 = FastChain((x, p) -> x.^3,
                  FastDense(2, 50, tanh),
                  FastDense(50, 2))
neural_ode_f(u,p,t) = dudt2(u,p)
pinit = initial_params(dudt2)
prob = ODEProblem(neural_ode_f, u0, tspan, pinit)

function predict_neuralode(p)
  tmp_prob = remake(prob,p=p)
  Array(solve(tmp_prob,Tsit5(),saveat=tsteps,sensealg=BacksolveAdjoint(autojacvec=ReverseDiffVJP(true))))
end

function loss_neuralode(p)
    pred = predict_neuralode(p)
    loss = sum(abs2, ode_data .- pred)
    return loss, pred
end

callback = function (p, l, pred; doplot = true)
  #display(l)
  # plot current prediction against data
  #plt = scatter(tsteps, ode_data[1,:], label = "data")
  #scatter!(plt, tsteps, pred[1,:], label = "prediction")
  #if doplot
  #  display(plot(plt))
  #end
  return false
end

@time result_neuralode = DiffEqFlux.sciml_train(loss_neuralode, pinit,
                                          ADAM(0.05), cb = callback,
                                          maxiters = 500)

#=
2.687161 seconds (17.79 M allocations: 1002.418 MiB, 7.41% gc time)

* Status: success

* Candidate solution
  Final objective value:     2.761669e-02

* Found with
  Algorithm:     ADAM

* Convergence measures
  |x - x'|               = NaN ≰ 0.0e+00
  |x - x'|/|x'|          = NaN ≰ 0.0e+00
  |f(x) - f(x')|         = NaN ≰ 0.0e+00
  |f(x) - f(x')|/|f(x')| = NaN ≰ 0.0e+00
  |g(x)|                 = NaN ≰ 0.0e+00

* Work counters
  Seconds run:   3  (vs limit Inf)
  Iterations:    500
  f(x) calls:    500
  ∇f(x) calls:   500
=#