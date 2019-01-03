using Plots
using Statistics
pyplot()
include("read_data.jl")
include("dmdc.jl")
include("helpers.jl")

# Set window ranges
dfile = "A_B_Uhat.h5"
dir = "prop_control/"
data_index = Colon()
w = 500
stp = 7
T = 4 # prediction window

# Fill the data
println("Loading data from files...")
data_dict = h5_to_dict(get_filename(1, dir))
q = length(data_dict["control_input"])
n = length(data_dict["sol_data"])
(data_index != Colon()) && (n = div(n,4))
Omega = Array{Float64,2}(undef, n+q, w-1)
Xp = Array{Float64,2}(undef, n, w-1)
fill_control_snapshots(Omega, Xp, dir, data_index)
X = Omega[1:n, :]
ind1, exact_sol1 = exact_sol_at_point(X)
println("done")

# load dynamics
f = h5open(dfile, "r")
A = read(f, "A")
B = read(f, "B")
U_hat = read(f, "U_hat")

# Compute the average running error
avg_err = average_prediction_error(A, B, U_hat, Omega, 1:stp:w-T-1, T)

# plot the results
p = plot(1:stp:w-T-1, log.(avg_err))
title!(string("Average Prediction Error \ndynamics: ", dfile, "\ndata: ", dir))
xlabel!("Iteration")
ylabel!("log(error)")
display(p)
savefig("A_from_prop4x_700_on_dmdc_30Dec_data")

