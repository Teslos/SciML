using DifferentialEquations
using Plots
using LinearAlgebra

function interaction(B, i, p)
    return 1 - dot(B, p[:A][i,:]) / p[:K][i]
end

function lotka_volterra!(dB, B, par, t)
    for i = 1:length(B)
        dB[i] = par[:r][i] * B[i] * interaction(B, i, par)
    end
end

A =[ 1 1.09 1.52 0;
     0 1 0.44 1.36;
     2.33 0 1 0.47;
     1.21 0.51 0.35 1]
r = [1, 0.72, 1.53, 1.27]
K = ones(4)

params = Dict(:A => A, :r => r, :K => K)
B0 = [0.1, 0.5, 0.05, 0.4]
tspan = (0.0, 1000.0)
prob = ODEProblem(lotka_volterra!, B0, tspan, params)
sol = solve(prob)
plot(sol, lw=3)
