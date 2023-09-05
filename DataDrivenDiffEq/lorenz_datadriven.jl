# Generating some data by solving a differential equation
using DataDrivenDiffEq
using DataDrivenSparse
using ModelingToolkit
using OrdinaryDiffEq
using Plots
using LinearAlgebra

# Define the system of differential equations for the Lorenz system
function lorenz(u, p, t)
    x, y, z = u
    ẋ = 10.0 * (y - x)
    ẏ = x * (28.0 - z) - y
    ż = x * y - (8.0 / 3.0) * z
    return [ẋ, ẏ, ż]
end

# Define the initial conditions
u0 = [1.0; 0.0; 0.0]
tspan = (0.0, 100.0)
dt = 0.1
prob =  ODEProblem(lorenz, u0, tspan)
sol = solve(prob, Tsit5(), dt=dt, saveat=dt)

# start the data-driven model
ddprob = DataDrivenProblem(sol)

@variables t x(t) y(t) z(t)
u = [x; y; z]
basis = Basis(polynomial_basis(u,5), u, iv = t)
opt  = STLSQ(exp10.(-5:0.1:-1))
ddsol = solve(ddprob, basis, opt, options = DataDrivenCommonOptions(digits = 1))

println(get_basis(ddsol))
