function kuramoto2d(N, T, t0, dt, s, K0, K1, sd_omega)
    s_sqrt_dt = sqrt(dt) * s
    omega = 1.0 .+ sd_omega * (randn(N, N) .- 0.5)
    phi = 2π * rand(N,N) .- π
    X = zeros(T, N,N)
    dphi = zeros(N,N)
    phi_up = zeros(N,N)
    phi_down = zeros(N,N)
    phi_left = zeros(N,N)
    phi_right = zeros(N,N)

    # coupling constant time course
    u = K0*ones(t0) 
    v = [ _ for _ in range(K0, stop=K1, length=T-t0) ]
    K = vcat(u,v)

    # iterate
    for t in 1:T+t0
                phi_up = circleshift(phi, [-1,0])
                phi_down = circleshift(phi, [1, 0])
                phi_left = circleshift(phi, [0, -1])
                phi_right = circleshift(phi, [0, 1])
                dphi = zeros(N,N)
                dphi += (sin.(phi_up .- phi) .+ sin.(phi_down .- phi) .+ sin.(phi_left .- phi) .+ sin.(phi_right .- phi))
                dphi = omega .+ K4[t] .* dphi
                phi += dt .* dphi + s_sqrt_dt * randn(N,N)
                (t > t0) && (X[t-t0,:,:] = angle.(exp.(im .* phi)))
            end
        println("\n")
        return X 
    end
