using OrdinaryDiffEq, SciMLSensitivity

function lotka_volterra(du,u,p,t)
    x, y = u
    α, β, δ, γ = p
    du[1] = dx = α*x - β*x*y
    du[2] = dy = -δ*y + γ*x*y
end

p = [1.5, 1.0, 3.0, 1.0]
prob = ODEProblem(lotka_volterra,[1.0;1.0],(0.0,10.0),p)
sol = solve(prob,Vern9(), abstol = 1e-10, reltol = 1e-10)

dg(out, u, p, t, i) = (out .= 1.0 .- u)
# if we had data we could do data[i] - u
ts = 0:0.5:10
res = adjoint_sensitivities(sol, Vern9(), dgdu_discrete = dg, t = ts, 
    abstol = 1e-10, reltol = 1e-10)
using ForwardDiff, Calculus, ReverseDiff, Tracker

function G(p)
    tmp_prob = remake(prob, u0 = convert.(eltype(p), prob.u0), p = p)
    sol = solve(tmp_prob, Vern9(), abstol = 1e-10, reltol = 1e-10,saveat=ts,
    sensealg = SensitivityADPassThrough())
    A = convert(Array,sol)
    sum((1 .- A).^2 ./ 2)
end
res2 = ForwardDiff.gradient(G, p)