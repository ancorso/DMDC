using Plots
pyplot()
include("cyl_data.jl")

last_snap = 1000
q,n = length.(read_snapshot(get_controlled_cyl_filename(1)))
n = convert(Int, n / 4)
Omega = Array{Float64,2}(n+q, last_snap-1)
Xp = Array{Float64,2}(n, last_snap-1)
fill_control_snapshots(Omega, Xp)

plot(1:last_snap-1, Omega[end, :])


## Compute the singular value decomposition of Omega
U, S, V = svd(Omega)
thresh = 0.1
rtil = length(S[S .> thresh])

# Truncate the matrices
Util = U[:,1:rtil]
Stil = S[1:rtil]
Stil = diagm(Stil)
Vtil = V[:, 1:rtil]

# Compute this efficient SVD of the output space Xp
U, S, V = svd(Xp)
r = length(S[S .> thresh])
Uhat = U[:,1:r]
Shat = S[1:r]
Shat = diagm(Shat)
Vbar = V[:,1:r]

U_1 = Util[1:n, :]
U_2 = Util[n+1:n+q, :]

approxA = Uhat' * Xp * Vtil * inv(Stil) * U_1' * Uhat
approxB = Uhat' * Xp * Vtil * inv(Stil) * U_2'

D, W = eig(approxA)
phi = (Xp * Vtil*inv(Stil)) * (U_1' * Uhat * W)


ppx = reshape(real(phi[:,116]), 256, 128)[:,:]
plot(1:256, 1:128, ppx', fill=true)

