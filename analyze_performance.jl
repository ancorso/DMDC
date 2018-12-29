# Loads in solution files and plots the cost function and control input
# as a function of time.
using HDF5
using Plots
pyplot()
include("read_data.jl")

fldr = "dmdc_online_new_opt/"
min_file = 1
max_file = 756


control = []
cost = []
for iteration in range(min_file, max_file)
    println("Loading file: ", iteration)
    dict = h5_to_dict(get_filename(iteration, fldr))
    push!(cost, dict["cost"])
    push!(control, dict["control_input"])
end


x = range(min_file, max_file)
p1 = plot(x, cost, title = "Cost function vs. iteration", xlabel="Iteration", ylabel = "Cost")
p2 = plot(x, control, title = "Control function vs. iteration", xlabel="Iteration", ylabel = "Control")

plot(p1, p2, size = (1200,400))
savefig("dmdc_online_w128_partialA0_new_opt")
