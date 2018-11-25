include("cyl_data.jl")
include("utils.jl")

# Read in the snapshots
last_snap = 6092
println("Reading in snapshots 1:", last_snap)
q,n = length.(read_snapshot(get_controlled_cyl_filename(1)))
n = convert(Int, n / 4)
Omega = Array{Float64,2}(undef, n+q, last_snap-1)
Xp = Array{Float64,2}(undef, n, last_snap-1)
l, n = size(Omega, 2), size(Xp, 1)
u,x = read_snapshot(get_controlled_cyl_filename(1))
fill_control_snapshots(Omega, Xp)

# Compute the dynamics
println("Computing dynamics...")
A, B, phi, D = DMDc(Omega, Xp)
println("A was ", size(A))


# save the corresponding dynamics
println("Saving dynamics to file...A...")
h5write("dynamics.h5", "A", A)
print("B...")
h5write("dynamics.h5", "B", B)
print("phi...")
h5write("dynamics.h5", "phi", abs.(phi))
println("D...")
h5write("dynamics.h5", "D", abs.(D))
