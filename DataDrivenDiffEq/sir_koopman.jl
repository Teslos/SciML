using OrdinaryDiffEq
using DataDrivenDiffEq
using DataDrivenDMD
using ModelingToolkit

function sir_ode(u,p,t)
    (s,i,r) = u
    (β,γ) = p
    ds = -β*s*i
    di = β*s*i - γ*i
    dr = γ*i
    [ds,di,dr]
end

p = [0.5,0.25]
u0 = [0.99, 0.01, 0.0]
tspan = (0.0, 40.0)
solver = ExplicitRK()
sir_prob = ODEProblem(sir_ode, u0, tspan, p)
sir_sol = solve(sir_prob, solver)
dd_prob = ContinuousDataDrivenProblem(sir_sol, sir_sol.t)
@parameters t
@variables (u(t))[1:3]
Ψ = Basis([u; u[1]*u[2]], u, independent_variable = t)
res_koopman = solve(dd_prob, Ψ, DMDPINV(), digits = 2)
sys_koopman = get_basis(res_koopman)
sys_parameters = get_parameter_map(sys_koopman)
