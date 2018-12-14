using HDF5
using Plots
pyplot()
include("cyl_data.jl")

fldr = "dfc_run_data/"
min_file = 1
max_file = 755

anim = @animate for iteration in range(min_file, max_file)
    println("Loading file: ", iteration)
    dict = h5_to_dict(get_filename(iteration, fldr))
    sol_data = dict["sol_data"][3,:,:]
    p = plot(1:256, 1:128, sol_data', fill=true)
end


gif(anim, string(rstrip(fldr, ['/']), ".gif"), fps=25)


