using OrdinaryDiffEq
using Plots

initial = [0.,0.1,0.5,0]
tspan = (0.,100.)

# Total energy of the system is consists
V(x,y) = 1//2 * (x^2 + y^2 + 2x^2*y - 2//3 * y^3)
E(x,y,dx,dy) = V(x,y)+1//2*(dx^2 + dy^2)

# define the function
function Henon_Heiles(du,u,p,t)
    x = u[1]
    y = u[2]
    dx = u[3]
    dy = u[4]

    du[1] = dx
    du[2] = dy
    du[3] = -x - 2x*y
    du[4] = y^2 - y - x^2

end

#Pass to solvers
prob = ODEProblem(Henon_Heiles, initial, tspan)
sol = solve(prob, Vern9(), abs_tol = 1e-16, rel_tol = 1e-16)

# Plot orbit
plot(sol, vars=(1,2), title="The orbit of the Henon Heiles system", xaxis="x", yaxis="y", leg=false)

plot(sol, vars=(1,3), title="Phase space for the Henon Heiles system", xaxis= "Position", yaxis="Velocity")
plot!(sol, vars=(2,4), leg=false)