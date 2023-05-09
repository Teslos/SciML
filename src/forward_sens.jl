using OrdinaryDiffEq, SciMLSensitivity

function lotka_volterra(du,u,p,t)
    x, y = u
    α, β, δ, γ = p
    du[1] = dx = α*x - β*x*y
    du[2] = dy = -δ*y + γ*x*y
end

p = [1.5, 1.0, 3.0, 1.0]
prob = ODEForwardSensitivityProblem(lotka_volterra,[1.0;1.0],(0.0,10.0),p)

sol = solve(prob,DP8())

x, dp = extract_local_sensitivities(sol)
da = dp[1]
x, dp = extract_local_sensitivities(sol,1)
using Plots
plot(sol.t,da',lw=3)
