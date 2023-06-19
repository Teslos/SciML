# Model to chaotic trains

function chaotictrains(du, u, params, t)
    l1 = params[1]
    l2 = params[2]
    m1 = params[3]
    m2 = params[4]
    g  = params[5]

    p1 = u[1]
    p2 = u[2]
    θ1 = u[3]
    θ2 = u[4]
    


    du .= [d_p1, d_p2, d_θ1, d_θ2]
    return nothing
end