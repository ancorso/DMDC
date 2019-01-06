include("read_data.jl")
include("dmdc.jl")
include("helpers.jl")

# Set window ranges
dir = "sin_omega/"
data_index = Colon()
ts = 1 # training starting location
te = 600 # training end
thresh = 1e-9 # threshold for dmdc

if length(ARGS) > 1
    println(ARGS)
    ts = parse(Int, ARGS[1])
    te = parse(Int, ARGS[2])
end
w = te - ts

# Load in the data from the files
println("Loading in data from files....")
data_dict = h5_to_dict(get_filename(1, dir))
q = length(data_dict["control_input"])
n = length(data_dict["sol_data"])
(data_index != Colon()) && (n = div(n,4))
Omega = Array{Float64,2}(undef, n+q, w)
Xp = Array{Float64,2}(undef, n, w)
# TODO:update this to fill in the appropriate snapshots
fill_control_snapshots(Omega, Xp, dir, data_index, ts, te)

X = Omega[1:n, :]
ind1, exact_sol1 = exact_sol_at_point(X)

println("done")

# perform the dynamic mode decomposition with control
# A, B, phi, W, U_hat= DMDc(Omega[:,1:tw], Xp[:,1:tw], 1e-6)
println("Performing dmdc...")
A, B, phi, W, U_hat = DMDc(Omega, Xp, thresh)
println("done")

# Save the results to the file
h5open("A_B_Uhat.h5", "w") do file
    write(file, "A", A)
    write(file, "B", B)
    write(file, "U_hat", U_hat)
end

# ppx = reshape(real(phi[:,1]), 256, 128)[:,:]
# plot(1:256, 1:128, ppx', fill=true)


