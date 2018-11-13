using Plots
pyplot()
include("cyl_data.jl")
include("utils.jl")

last_snap = 1000
q,n = length.(read_snapshot(get_controlled_cyl_filename(1)))
n = convert(Int, n / 4)
Omega = Array{Float64,2}(undef, n+q, last_snap-1)
Xp = Array{Float64,2}(undef, n, last_snap-1)
l, n = size(Omega, 2), size(Xp, 1)
u,x = read_snapshot(get_controlled_cyl_filename(1))
fill_control_snapshots(Omega, Xp)

plot(1:last_snap-1, Omega[end, :])


A, B, phi, D = DMDc(Omega, Xp)


ppx = reshape(real(phi[:,116]), 256, 128)[:,:]
plot(1:256, 1:128, ppx', fill=true)

