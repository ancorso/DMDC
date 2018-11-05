using Plots
using HDF5

const rho = 1
const ux = 2
const uy = 3
const p = 4

h5read()
f = h5open("data/rot_cyl_re50/sol_data_0001.h5", "r")
names(f)
u = read(f, "control_input")
x = read(f, "sol_data")



plot(1:256, 1:128, x[4,:,:]')

