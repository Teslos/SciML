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

g(u,p,t) = (sum(u) .^ 2) ./ 2
function dg(out, u, p, t)
    out[1] = u[1] + u[2]
    out[2] = u[1] + u[2]
end

prob = ODEProblem(lotka_volterra,[1.0;1.0],(0.0,10.0),p)
sol = solve(prob,DP8())
res = adjoint_sensitivities(sol, Vern9(), g = g, dgdu_continuous = dg, sensealg=BacksolveAdjoint())

using QuadGK, ForwardDiff, Calculus
function G(p)
    tmp_prob = remake(prob, u0 = convert.(eltype(p), prob.u0), p = p)
    sol = solve(tmp_prob, Vern9(), abstol = 1e-10, reltol = 1e-10)
    res, err = quadgk((t)->(sum(sol(t)).^2) ./ 2, 0, 10, atol = 1e-14, rtol = 1e-14)
    res
end
res2 = ForwardDiff.gradient(G, p)
res3 = Calculus.gradient(G, p)