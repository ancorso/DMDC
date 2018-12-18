# Compute the dynamics matrices and eigenvectors from solution data outputted
# directly from PyFR. this is used to contruct the initial A,B,and U, or those
# used for the offline DMDc

include("read_data.jl")
include("dmdc.jl")

last_snapshot = 600
dir = "data/re50_17Dec2018/"

# Read the number of snapshots from the command line
if length(ARGS) > 0
    last_snapshot = convert(Int, ARGS[1])
end
if length(ARGS) > 1
    dir = ARGS[2]
end

# Read in the snapshots
println("Reading in snapshots 1:", last_snapshot)
q,n = length.(read_snapshot(get_filename(1,dir)))

Omega = Array{Float64,2}(undef, n+q, last_snapshot-1)
Xp = Array{Float64,2}(undef, n, last_snapshot-1)
l, n = size(Omega, 2), size(Xp, 1)
fill_control_snapshots(Omega, Xp, dir)

# Compute the dynamics
println("Computing dynamics...")
A, B, phi, D, U_hat = DMDc(Omega, Xp, 1e-6)
U_hat = convert(Array{Float64, 2}, U_hat')
println("A was ", size(A))
println("Uhat is: ", typeof(U_hat), " size: ", size(U_hat))

# save the corresponding dynamics
println("Saving dynamics to file...")
h5open("A_B_Uhat.h5", "w") do file
    write(file, "A", A)
    write(file, "B", B)
    write(file, "U_hat", U_hat)
end

