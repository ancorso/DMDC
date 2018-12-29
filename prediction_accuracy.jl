using Plots
using Statistics
pyplot()
include("read_data.jl")
include("dmdc.jl")
include("helpers.jl")

#note that A_B_Uhat_good.hdf5 has the params (ts = 1, tw = 300, ps = 100, pw = 30, w = 999)

# Set window ranges
ts = 1 # training starting location
tw = 900 # training window
ps = 501
# ps = ts + tw
pw = 30 # prediction window
w = 999 # all the data we need
thresh = 1e-9 # threshold for dmdc
recompute_dmd = false
refill_arrays = true

if refill_arrays
    println("Loading in data from files....")
    # Read in the data
    dir = "data/re50_17Dec2018/"
    data_index = 3
    q,n = length.(read_snapshot(get_filename(1, dir)))
    # n, = length.(read_snapshot(get_filename(1, dir)))

    n = div(n,4)
    Omega = Array{Float64,2}(undef, n+q, w-1)

    # X = Array{Float64,2}(undef, n, w-1)
    Xp = Array{Float64,2}(undef, n, w-1)
    fill_control_snapshots(Omega, Xp, dir, data_index)
    # fill_static_snapshots(X, Xp, dir)

    X = Omega[1:n, :]
    ind1, exact_sol1 = exact_sol_at_point(X)

    println("done")
end



# perform the dynamic mode decomposition with control
# A, B, phi, W, U_hat= DMDc(Omega[:,1:tw], Xp[:,1:tw], 1e-6)
if recompute_dmd
    println("Performing dmdc...")
    # A, phi, D, W, U_hat = DMD(X[:,1:tw], Xp[:,1:tw], thresh)
    A, B, phi, W, U_hat = DMDc(Omega[:,ts:ts+tw-1], Xp[:,ts:ts+tw-1], thresh)
    println("done")

    println("renormalizing the dynamics matrix...")

    D, W = eigen(A)
    A = real.(W * diagm(0=> cap.(D)) * inv(W))
    println("done")
end


# Find the point of largest total error
println("finding max error point...")
x = X[:,ps-1]
spatial_errors = zeros(size(x))
for i=ps:ps + pw - 1
    xtil = U_hat' * x
    # Xtil_p1 = A * xtil
    Xtil_p1 = A * xtil + B * Omega[n+1:n+q, i]
    global x = U_hat * Xtil_p1

    global spatial_errors += abs.(X[:, i] - x)
end
max_val, ind2 = findmax(spatial_errors)
println("done")



# Approximate the solution
x = X[:,ps-1]
errl = []
approx_sol1 = []
approx_sol2 = []
exact_sol2 = []
for i=ps:ps + pw - 1
    println("predicting: ", i)
    xtil = U_hat' * x
    Xtil_p1 = A * xtil + B * Omega[n+1:n+q, i]
    global x = U_hat * Xtil_p1

    push!(errl, norm(X[:, i] - x))
    push!(approx_sol1, x[ind1])
    push!(approx_sol2, x[ind2])
    push!(exact_sol2, X[ind2,i])
end

println("Average error: ", mean(errl))
p1 = plot(exact_sol1[ps:ps+pw-1], label="exact solution")
plot!(approx_sol1, label="approximate solution")
plot!(Omega[n+q,ps:ps+pw-1], label="control input")
title!("Max Dynamics Point")

p2 = plot(exact_sol2, label="exact solution")
plot!(approx_sol2, label="approximate solution")
title!("Max Error Point")


p = plot(p1, p2)

display(p)
plot(X[13999, 1:998])
plot!(Omega[n+1, 1:998])

h5open("A_B_Uhat.h5", "w") do file
    write(file, "A", A)
    write(file, "B", B)
    write(file, "U_hat", U_hat)
end

println("max dynamics index: ", to_ijk(ind1, [256, 128] ))
println("max err index: ", to_ijk(ind2, [256, 128] ))

D, W = eigen(A)
sort(abs.(D))

plot(1:256, 1:128, reshape(X[:,1], 256, 128)', fill=false)


