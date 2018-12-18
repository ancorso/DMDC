# Compute the dynamics matrices and eigenvectors from solution data outputted
# directly from PyFR. this is used to contruct the initial A,B,and U, or those
# used for the offline DMDc

include("read_data.jl")
include("dmdc.jl")

# Read in the snapshots
f = h5open("X_u.h5", "r")
X = read(f, "X")
u = read(f, "u")
T = size(X,4)
X = reshape(X, (:, T))

# concatenate the data to form the input and output matrices
size(X[:, 1:T-1])
Omega = vcat(X[:, 1:T-1], u[:,1:end])
Xp = X[:, 2:T]

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

