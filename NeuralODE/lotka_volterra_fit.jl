using DifferentialEquations, Optimization, OptimizationPolyalgorithms, SciMLSensitivity
# lotka-volterra function
function lotka_volterra!(du, u, p, t)
    x, y = u
    α,β, δ, γ = p
    du[1] = α * x - β * x * y
    du[2] = -δ * y + γ * x * y
end

# initial condition 
u0 = [1.0, 1.0]
# Simulation interval 
tspan = (0.0, 10.0)

# LV equation parameters p = [α, β, δ, γ]
p = [1.5, 1.0, 3.0, 1.0]

# Setup the ODE problem, then solve
prob = ODEProblem(lotka_volterra!, u0, tspan, p)
datasol = solve(prob, saveat = 1)

function loss(newp)
    newprob = remake(prob, p = newp)
    sol = solve(newprob, saveat = 1)
    loss = sum(abs2, sol .- datasol)
    return loss, sol
end

callback = function(p,l,sol)
    display(l)
    plt = plot(sol, ylim = (0,6), label="Current prediction")
    scatter!(plt, datasol, label="Data")
    display(plt)
    # tell optimization solve to halt optimization. if return true, 
    # then optimization stops.
    return false
end

adtype = Optimization.AutoForwardDiff()
optf  = Optimization.OptimizationFunction((x,p) -> loss(x), adtype)
pguess = [1.0, 1.2, 2.5, 1.2]
optprob = Optimization.OptimizationProblem(optf, pguess)

result_ode = Optimization.solve(optprob, PolyOpt(),
                                callback = callback,
                                maxiters = 200)

