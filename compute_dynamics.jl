include("cyl_data.jl")
include("utils.jl")

# Read in the snapshots
# last_snap = 6092
last_snap = 100
println("Reading in snapshots 1:", last_snap)
q,n = length.(read_snapshot(get_controlled_cyl_filename(1)))

Omega = Array{Float64,2}(undef, n+q, last_snap-1)
Xp = Array{Float64,2}(undef, n, last_snap-1)
l, n = size(Omega, 2), size(Xp, 1)
u,x = read_snapshot(get_controlled_cyl_filename(1))
fill_control_snapshots(Omega, Xp)

# Compute the dynamics
println("Computing dynamics...")
A, B, phi, D, U_hat = DMDc(Omega, Xp)
U_hat = convert(Array{Float64, 2}, U_hat')
println("A was ", size(A))
println("Uhat is: ", typeof(U_hat), " size: ", size(U_hat))

# save the corresponding dynamics
println("Saving dynamics to file...")
h5open("dynamics.h5", "w") do file
    write(file, "A", A)
    write(file, "B", B)
    write(file, "U_hat", U_hat)
end

