using DifferentialEquations
using Plots
p = (1.0, 1.0, 1.0, 10.0, 0.001, 100.0) # a, α, ubar, β, D1, D2 
N = 100
Ax = Tridiagonal([1.0 for i in 1:N-1], [-2.0 for i in 1:N], [1.0 for i in 1:N-1])
Ay = copy(Ax)
Ax[2,1] = 2.0
Ax[end-1,end] = 2.0
Ay[1,2] = 2.0
Ay[end,end-1] = 2.0

function basic_version!(dr,r,p,t)
    a, α, ubar, β, D1, D2 = p
    u = r[:,:,1]
    v = r[:,:,2]
    Du = D1*Ay*u + D1*u*Ax
    Dv = D2*Ay*v + D2*v*Ax
    dr[:,:,1] = Du .+ a*u*u./v .+ ubar -α*u
    dr[:,:,2] = Dv .+ a*u*u .- β*v
end

a, α, ubar, β, D1, D2 = p
uss = (ubar+β)/a
vss = (a*uss^2)/β
r0 = zeros(N,N,2)
r0[:,:,1] .= uss+0.1*rand()
r0[:,:,2] .= vss
prob = ODEProblem(basic_version!, r0, (0.0, 0.1), p)
solve(prob, Tsit5(), saveat=0.01)