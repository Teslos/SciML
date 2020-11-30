using DiffEqFlux, OrdinaryDiffEq, Flux, Optim, Plots
u0 = Float32[2.0;0.0]
tspan = (0.0f0, 1.5f0)
datasize = 30
tsteps = range(tspan[1],tspan[2],length=datasize)
function trueODEfunc(du,u,p,t)
    true_A = [-0.1 2.0; -2.0 -0.1]
    du .= ((u.^3)'true_A)'
end

prob_trueode = ODEProblem(trueODEfunc,u0,tspan)
ode_data = Array(solve(prob_trueode,Tsit5(),saveat=tsteps))

dudt2 = FastChain((x,p) -> x.^3,
                FastDense(2,50,tanh),
                FastDense(50,2))
neural_ode_f(u,p,t) = dudt2(u,p)

pinit = initial_params(dudt2)
prob = ODEProblem(neural_ode_f,u0,(0.0f0,1.5f0),pinit)
sol = solve(prob, saveat=tsteps)

plot(sol)
scatter!(tsteps,ode_data')

function loss(p; salg=InterpolatingAdjoint())
    tmp_prob = remake(prob,p=p)
    tmp_sol = solve(tmp_prob, Tsit5(), saveat=tsteps, sensealg=salg)
    sum(abs2, Array(tmp_sol)-ode_data)
end

function neuralode_callback(p,l)
    @show l
    tmp_prob = remake(prob,p=p)
    tmp_sol = solve(tmp_prob, Tsit5(), saveat=tsteps)
    fig = plot(tmp_sol)
    scatter!(fig,ode_data')
    display(fig)
    false
end

@time res = DiffEqFlux.sciml_train(p->loss(p;salg=BacksolveAdjoint(autojacvec=ReverseDiffVJP(true))), pinit, ADAM(0.05),
    maxiters=100,cb=neuralode_callback)

@time DiffEqFlux.sciml_train(loss,res.minimizer,BFGS(initial_stepnorm=0.01),
    cb=neuralode_callback)

@time DiffEqFlux.sciml_train(p->loss(p;salg=BacksolveAdjoint(autojacvec=ReverseDiffVJP(true))),
        res.minimizer, BFGS(initial_stepnorm=0.01),
        maxiters=100,cb=neuralode_callback)
