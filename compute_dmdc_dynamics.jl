# Compute the dynamics matrices and eigenvectors from solution data outputted
# directly from PyFR.
# This is used to contruct the initial A,B,and U, or those used for offline DMDc
# usage:
#           1) julia compute_dmdc_dynamics.jl matrix filename
#           2) julia compute_dmdc_dynamics.jl files dir first_iter last_iter
include("dmdc.jl")

export_type = ARGS[1]

# Non-input parameters
data_index = Colon()
thresh = 0.99 # Fraction of energy retained in SVD
output_filename = "A_B_Uhat.h5"
num_modes_override = nothing

# Load in the data by whatever means specified
print("Loading data...")
if export_type == "matrix"
    filename = ARGS[2]
    f = h5open(filename, "r")
    X, u = read(f, "X"), read(f, "u")
    close(f)
    X = reshape(X[data_index, :, :, :], (:, size(X, 4)))
    global Ω, Xp = vcat(X[:, 1:end-1], u[:,1:end]), X[:, 2:end]
elseif export_type == "files"
    include("load_data.jl")
    dir, tstart, tend = ARGS[2], parse(Int, ARGS[3]), parse(Int, ARGS[4])
    if length(ARGS) == 5
        print("Setting mode count override to ", ARGS[5], "...")
        num_modes_override = parse(Int, ARGS[5])
    end
    global Ω, Xp = load_data(dir, tstart:tend, data_index = data_index)
else
    @error "Invalid PyFR export type"
end

# Compute the dynamics
print("Computing DMDc...")
A, B, phi, D, U_hat = DMDc(Ω, Xp, thresh, num_modes_override = num_modes_override)

# save the corresponding dynamics
print("Saving file to \"", output_filename, "\"...")
h5open(output_filename, "w") do file
    write(file, "A", A)
    write(file, "B", B)
    write(file, "U_hat", U_hat)
end
println("done!")