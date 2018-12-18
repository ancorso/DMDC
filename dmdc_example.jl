# Example of a dynamic mode decomposition with control

using Plots
pyplot()
include("read_data.jl")
include("dmdc.jl")

last_snap = 1000
q,n = length.(read_snapshot(get_filename(1, "data/re50_17Dec2018/")))
n = convert(Int, n/4)
Omega = Array{Float64,2}(undef, n+q, last_snap-1)
Xp = Array{Float64,2}(undef, n, last_snap-1)
l, n = size(Omega, 2), size(Xp, 1)
u,x = read_snapshot(get_filename(1, "data/re50_17Dec2018/"))
fill_control_snapshots(Omega, Xp)

plot(1:last_snap-1, Omega[end, :])


A, B, phi, W, U_hat = DMDc(Omega, Xp)

ppx = reshape(real(phi[:,end-4]), 256, 128)
plot(1:256, 1:128, ppx', fill=false)

savefig("m1")

