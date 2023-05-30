using Flux, DiffEqFlux, DifferentialEquations, Statistics, Plots, ReverseDiff
t = range(0.0f0, 1.0f0, length = 1024)
π_32 = Float32(π)

qₜ = reshape(sin.(2π_32 * t), 1, :)
pₜ = reshape(cos.(2π_32 * t), 1, :)
dqdt = 2π_32 .* pₜ
dpdt = -2π_32 .* qₜ

data = cat(qₜ, pₜ, dims=1)
target = cat(dqdt, dpdt, dims = 1)
dataloader = Flux.Data.DataLoader((data, target); batchsize=256, shuffle = true)

hnn = HamiltonianNN(
    Flux.Chain(Flux.Dense(2,64,relu), Flux.Dense(64,1))
)

p = hnn.p 
opt = ADAM(0.01)

# define the loss function
loss(x, y, p) = mean((hnn(x,p) .- y) .^ 2)
mycallback() = println("Loss Neural Hamiltonian  DE = $(loss(data, target, p))")

epochs = 500

for epoch in 1:epochs
    for (x,y) in dataloader
        gs = ReverseDiff.gradient(p -> loss(x,y,p), p)
        Flux.Optimise.update!(opt, p, gs)
    end
    if epoch % 100 == 1
        mycallback()
    end
end
mycallback()

model = NeuralHamiltonianDE(
    hnn, (0.0f0, 1.0f0),
    Tsit5(), save_everystep = false,
    save_start = true, saveat = t
)

pred = Array(model(data[:,1]))
plot(data[1,:],data[2,:], lw=4, label="Original")
plot!(pred[1,:],pred[2,:], lw=4, label="Predicted")
xlabel!("Position (q)")
ylabel!("Momentum (p)")



