using HDF5
using Plots
pyplot()
include("read_data.jl")

fldr = "dmdc_30Dec_A_from_prop/"
min_file = 1
max_file = 754

anim = @animate for iteration in range(min_file, max_file)
    println("Loading file: ", iteration)
    dict = h5_to_dict(get_filename(iteration, fldr))
    sol_data = dict["sol_data"][3,:,:]
    p = plot(1:256, 1:128, sol_data')
    title!(string("Y-Vel at Iteration: ", iteration))
end


gif(anim, string(rstrip(fldr, ['/']), ".gif"), fps=25)


