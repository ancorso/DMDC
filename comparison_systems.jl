include("dmdc.jl")
using Plots
pyplot()
using Statistics

# setup the sample system
T, n, q = 50, 10, 1

max_iter = 1000
iter  = 0
while iter < max_iter
    global iter += 1
    global A = rand(n,n) .- 0.5
    global D, W = eigen(A)

    if all(abs.(D) .<= 1) break end
end
D

A = diagm(0=>[.1, .01, .01, .01, .01, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6 ])
B = rand(n, q) .- 0.5

u = rand(T)
x0 = rand(n)

# Fill the data snapshots
Omega = Array{Float64, 2}(undef, n+q, T-1)
Xp = Array{Float64, 2}(undef, n, T-1)
for i=1:T-1
    Omega[1:n, i] = x0
    Omega[n+1:n+q, i] .= u[i]
    global x1 = A*x0 + B*u[i]
    Xp[1:n, i] = x1
    global x0 = x1
end

# Compute the dmdc
At, Bt, phi, W, U_hat = DMDc(Omega, Xp, 1e-3)
At

x = Omega[1:n,1]
err_l = []
for i=2:T-1
    println("predicting: ", i)
    xtil = U_hat' * x
    Xtil_p1 = At * xtil + Bt*Omega[n+1:n+q, i-1]
    global x = U_hat * Xtil_p1

    # Get the stepwise error
    err = abs.(Omega[1:n, i] - x)
    push!(err_l, mean(err))
end

err_l
plot(err_l, label="error")

