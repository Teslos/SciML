using DataDrivenDiffEq
using DataDrivenSparse
using ModelingToolkit
using OrdinaryDiffEq

using LinearAlgebra
using Plots

# Create a test problem
function lorentz(u,p,t)
    x, y, z = u
    dx_dt = 10*(y-x)
    dy_dt = x*(28.0-z)-y
    dz_dt = x*y - (8/3)*z
    return [dx_dt,dy_dt,dz_dt]
end

u0 = [8.0; 7.0; 27.0]
p = [10.0; -10.0; 28.0; -1.0; 1.0; -8/3]
tspan = (0.0,100.0)
dt = 0.001
prob = ODEProblem(lorentz, u0, tspan)
solution = solve(prob, Tsit5(), saveat = dt)

# @variables x y z
# X = Array(solution)
# DX = similar(X)
# for (i, xi) in enumerate(eachcol(X))
#     DX[:,i] = lorentz(xi, [], 0.0)
# end

# u = Operation[x; y; z]
# polys = Operation[]
# for i in 0:4
#     for j in 0:i
#         for k in 0:j
#             push!(polys, u[1]^i * u[2]^j * u[3]^k)
#             push!(polys, u[2]^i * u[3]^j * u[1]^k)
#             push!(polys, u[3]^i * u[1]^j * u[2]^k)
#         end
#     end
# end


# opt = STRRidge(0.1)
# psi = SInDy(X,DX,basis, maxiter=100, opt=opt, normalize=true)
# print_equations(psi)
# get_error(psi)

### start automatic discovery
ddprob = DataDrivenProblem(solution)
@variables t x(t) y(t) z(t)
u = [x;y;z]
basis = Basis(polynomial_basis(u,5), u, iv=t)
opt = STLSQ(exp10.(-5:0.1:-1))
ddsol = solve(ddprob, basis, opt, options=DataDrivenCommonOptions(digits=1))
println(get_basis(ddsol))
