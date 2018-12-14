using Plots
using LinearAlgebra

pyplot()
include("cyl_data.jl")
include("utils.jl")

last_snap = 999
n = convert(Int, length(read_snapshot(get_static_cyl_filename(999))[1]))
X = Array{Float64,2}(undef, n, last_snap)
Xp = Array{Float64,2}(undef, n, last_snap)
fill_static_snapshots(X, Xp)

# Animation of data
# anim = @animate for i=1:last_snap-1
#     println(string("saving iteration ", i))
    ppx = reshape(X[:,i], 4, 256, 128)[1,:,:]
#     plot(1:256, 1:128, ppx', fill=true)
# end
# gif(anim, "./static_cyl.gif", fps = 15)

Atil, phi, D = DMD(X, Xp)
size(X, 1)

pts = [20, 50, 100, 1000, 15000]
Aerr = []
Derr = []
phierr = []
for p in pts
    Atilc, phic, Dc = compressed_DMD(X, Xp, p)
    push!(Aerr, sum(abs.(Atil - Atilc)) / length(Atil))
    push!(Derr, sum(abs.(D - Dc)) / size(D,1))
    push!(phierr, sum(abs.(phi - phic)) / length(phi))
end
Aerr

plot(pts, Aerr, label="A", yscale=:log10)
plot!(pts, Derr, label="D")
plot!(pts, phierr, label="phi")
title!("Compressed sensing error")
xlabel!("Number of sample points")
ylabel!("Average error")

ppx = reshape(real(phi[:,1]), 256, 128)[:,:]
plot(1:256, 1:128, ppx', fill=true)

println("hello")
# Time dynamics of data
# omega = log.(D)

# x1 = X[:,1]
# b = phi \ x1
# mm1 = size(X,2)
# time_dynamics = zeros(Complex, rtil, mm1)
# for iter = 1:mm1
#     time_dynamics[:,iter] = b.*exp.(omega*iter)
# end

# Xp = phi * time_dynamics

