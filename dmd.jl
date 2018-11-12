using Plots
pyplot()
include("cyl_data.jl")

last_snap = 999
n = convert(Int, length(read_snapshot(get_static_cyl_filename(999))[1]) / 4)
X = Array{Float64,2}(n, last_snap)
Xp = Array{Float64,2}(n, last_snap)
fill_static_snapshots(X, Xp)

# anim = @animate for i=1:last_snap-1
#     println(string("saving iteration ", i))
#     ppx = reshape(X[:,i], 4, 256, 128)[1,:,:]
#     plot(1:256, 1:128, ppx', fill=true)
# end
# gif(anim, "./static_cyl.gif", fps = 15)

# take the singular value decomposition
U, S, V = svd(X)
thresh = .1
rtil = length(S[S .> thresh])

# Truncate the matrices
Util = U[:,1:rtil]
Stil = diagm(S[1:rtil])
Vtil = V[:, 1:rtil]

sinv = inv(Stil)
Atil = conj(Util') * Xp * Vtil * sinv
D, W = eig(Atil)
phi = Xp * Vtil * sinv * W


omega = log.(D)

x1 = X[:,1]
b = phi \ x1
mm1 = size(X,2)
time_dynamics = zeros(Complex, rtil, mm1)
for iter = 1:mm1
    time_dynamics[:,iter] = b.*exp.(omega*iter)
end

Xp = phi * time_dynamics


ppx = reshape(real(phi[:,8]), 256, 128)[:,:]
plot(1:256, 1:128, ppx', fill=true)

