using Zygote, ForwardDiff

function s(x)
    t = 0.0
    sign = -1.0
    for i in 1:19
        if isodd(i)
            @show i
            newterm = x^i/factorial(i)
            abs(newterm) < 1e-8 && return t
            sign = -sign
            t += sign*newterm
        end
    end
    return t
end
ForwardDiff.derivative(s,1.0)
Zygote.gradient(s,1.0)

# definition of composition of the functions
# function J(f ∘ g)(x)
#     a,da = J(f)(x)
#     b,db = J(g)(x)
#     b,z -> da(db(z))
# end


f(x) = x^2 + 3x +1
Zygote.gradient(f,1/3)
