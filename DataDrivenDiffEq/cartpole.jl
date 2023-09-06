# cartpole example

using DataDrivenDiffEq
using ModelingToolkit
using OrdinaryDiffEq
using LinearAlgebra
using DataDrivenSparse

using Plots
gr()

# define the system
function cart_pole(u, p, t)
    du = similar(u)
    F = -0.2 + 0.5 * sin( 6 * t) # control input
    