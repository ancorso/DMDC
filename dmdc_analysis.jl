include("load_data.jl")
using Statistics
using LinearAlgebra
using Plots; gr()

# See what the percent difference that should be expected between frames
function expected_errs(X, n)
    w = size(X,2)
    max_avg_diff = 0
    min_avg_diff = 1e100
    for i=1:1000
        i1, i2 = rand(1:w-1), rand(1:w-1)
        if i1 == i2 break end

        err = norm(X[1:n, i1] - Omega[1:n, i2])
        if err < min_avg_diff
            min_avg_diff = err
        end
        if err > max_avg_diff
            max_avg_diff = err
        end
    end
    # max_avg_diff # ~ 14.7
    # min_avg_diff # ~ 1.2
    min_avg_diff, max_avg_diff
end

# Reads cost function and control input and plots both next to each other.
# dir - specifies the directory where the solution files are stored.
# iterations - specifies which iterations should be read in
# output_img_name - The name of the output image
function plot_suppression_performance(dir, iterations, output_img_name)
    # Load in cost and control input at desired iterations
    control, cost = get_scalar_sequences(dir, iterations, ["control_input", "cost"])

    # Plot the performance
    p1 = plot(iterations, cost, title = "Cost function vs. Iteration", xlabel="Iteration", ylabel = "Cost")
    p2 = plot(iterations, control, title = "Control function vs. Iteration", xlabel="Iteration", ylabel = "Control")
    plot(p1, p2, size = (1200,400))
    savefig(output_img_name)
end

# A, B and Uhat are the dmdc params
# Ω is the true data
# starting_points are the iterations that the prediction should start at
# T is the prediction window
function average_prediction_error(A, B, U_hat, Ω, starting_points, T)
    n = size(U_hat, 1)
    q = size(Ω,1) - n
    average_err = Float64[]
    for s in starting_points
        println("Predicting from: ", s)
        x = Ω[1:n,s]
        errl = Float64[]
        for i=s+1:s+T
            xtil = U_hat' * x
            Xtil_p1 = A * xtil + B * Ω[n+1:n+q, i-1]
            x = U_hat * Xtil_p1

            push!(errl, norm(Ω[1:n, i] - x))
        end
        push!(average_err, mean(errl))
    end
    average_err
end

# Plots the prediction accuracy over a specified range
function plot_prediction_accuracy(dynamics_file, dir, starting_points, T, output_img_name; read_control = true, data_index = Colon())
    # Load the dynamics
    print("Loading Dynamics...")
    f = h5open(dynamics_file, "r")
    A, B, U_hat = read(f, "A"), read(f, "B"), read(f, "U_hat")
    close(f)

    # Load the comparison data
    print("Loading Comparison Data...")
    Ω = load_data(dir, 1:maximum(starting_points) + T, read_control = read_control, data_index = data_index, get_next_frame = false)

    plot_prediction_accuracy(A, B, U_hat, Ω, starting_points, T, output_img_name)
end

# Plots the prediction accuracy over a specified range
function plot_prediction_accuracy(A, B, U_hat, Ω, starting_points, T, output_img_name)
    # Compute the average running error
    print("Running preditions:")
    avg_err = average_prediction_error(A, B, U_hat, Ω, starting_points, T)

    # Plot the results
    print("Plotting...")
    p = plot(starting_points, log.(avg_err), xlabel = "Iteration", ylabel="log(Average Error)", title = string("Average Prediction Error"))
    savefig(output_img_name)
    println("done!")
end

function make_side_by_side_vids()
    # fldr = "dmdc_30Dec_A_from_prop/"
    # min_file = 1
    # max_file = 754
    #
    # anim = @animate for iteration in range(min_file, max_file)
    #     println("Loading file: ", iteration)
    #     dict = h5_to_dict(get_filename(iteration, fldr))
    #     sol_data = dict["sol_data"][3,:,:]
    #     p = plot(1:256, 1:128, sol_data')
    #     title!(string("Y-Vel at Iteration: ", iteration))
    # end
    #
    #
    # gif(anim, string(rstrip(fldr, ['/']), ".gif"), fps=25)

end




