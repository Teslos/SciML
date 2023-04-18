using DifferentialEquations, Optimization, OptimizationPolyalgorithms, SciMLSensitivity, Zygote, Plots

function lotka_volterra!(du,u,p,t)
    x, y = u
    α, β, δ, γ = p
    du[1] = dx = α*x - β*x*y
    du[2] = dy = -δ*y + γ*x*y
end

# Initial parameters
u0 = [1.0,1.0]

# Simulation interval and intermediary points
tspan = (0.0,10.0)
tsteps = 0.0:0.01:10.0

# LV equation parameters p = [α, β, δ, γ]
p = [1.5, 1.0, 3.0, 1.0]


# Solve the ODE
prob = ODEProblem(lotka_volterra!,u0,tspan,p)
sol = solve(prob,Tsit5(),saveat=tsteps)

# Plot the solution
plot(sol,vars=(1,2),label=["x" "y"],title="Lotka-Volterra",lw=2)

# Define the cost function
function loss(p)
    sol = solve(prob,Tsit5(),saveat=tsteps,p=p)
    sum(abs2, sol .- 1)
    return loss, sol
end

callback = function(p, l, pred)
    display(l)
    plt = plot(pred, ylim = (0,6))
    display(plt)
    # Tell Optimization.solve to not halt the optimization.
    # If return true, Optimization.solve will halt the optimization.

    return false
end

# Define the optimization problem
adtype = Optimization.AutoZygote()
optf = Optimization.OptimizationFunction((x,p)->loss(x),adtype)
optprob = Optimization.OptimizationProblem(optf,p)

result_ode = Optimization.solve(optprob, PolyOpt(), callback=callback, maxiters=100)

# Plot the solution
remade_sol = solve(remake(sol,p=result_ode.u),Tsit5(),saveat=tsteps)
plot(remade_sol,vars=(1,2),label=["x" "y"],title="Lotka-Volterra",lw=2)