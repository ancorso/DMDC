using Plots
pyplot()

using HDF5

# Construct the h5 filename at the specified iteration
function get_filename(iteration)
    dir = "data/rot_cyl_re50/"
    base = "sol_data_"
    num = lpad(iteration, 4, 0)
    ext = ".h5"
    return string(dir, base, num, ext)
end

# Read the h5 file at the specified iteration and return the solution data and the control input
function get_snapshot(iteration)
    # Note that in each hdf5 file, there is density, ux, uy, energy
    f = h5open(get_filename(iteration), "r")
    names(f)
    u = read(f, "control_input")
    x = read(f, "sol_data")
    x = x[1,:,:]
    x,u
end

function get_snapshot_size(iteration = 1)
    x, u = get_snapshot(iteration)
    return length(x), length(u)
end

# Construct the solution data and control matrices
function fill_snapshots(Omega, Xp, first=1)
    l, n = size(Omega, 2)+1, size(Xp, 1)
    last = first+l-1
    for i=first:last
        (i % 100 == 0) && println(string("Adding snapshot ", i))
        x,u = get_snapshot(i)
        (i < last) && (Omega[1:n,i] = x[:])
        (i < last) && (Omega[n+1:end,i] = u)
        (i > 2) && (Xp[:,i-1] = x[:])
    end
    Omega, Xp
end

last_snap = 1000
n, q = get_snapshot_size()
Omega = Array{Float64,2}(n+q, last_snap-1)
Xp = Array{Float64,2}(n, last_snap-1)
fill_snapshots(Omega, Xp)

plot(1:last_snap-1, Omega[end, :])



## Compute the singular value decomposition of Omega
U, S, V = svd(Omega)
thresh = 0
S
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
x1 = abs.(reshape(phi[:,3], (1, 256, 128)))
plot(1:256, 1:128, x1[1,:,:]')


