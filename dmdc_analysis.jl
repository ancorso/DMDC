include("load_data.jl")
include("dmdc.jl")
using Statistics
using LinearAlgebra
using Plots; gr()

# Load the dynamics file
function load_dynamics(file)
    f = h5open(file, "r")
    A, B, U_hat = read(f, "A"), read(f, "B"), read(f, "U_hat")
    close(f)
    A,B,U_hat
end

function plot_modes(dynamics_file, output_img, num_modes = nothing)
    _,_,U_hat = load_dynamics(dynamics_file)
    r = (num_modes == nothing) ? size(U_hat, 2) : num_modes
    dofs = 4
    plots = []
    for dof = 1:dofs
        for m = 1:r
            mode = reshape(real(U_hat[:,m]), 4, 256, 128)
            push!(plots, plot(1:256, 1:128, mode[dof,:,:]', title = string("Mode: ", m, " dof: ", dof), size=(600,400)))
        end
    end
    plot(plots..., layout = (dofs, r), size = (600*r,400*dofs))
    savefig(output_img)
end

function plot_B(dynamics_file, output_img)
    _,B,U_hat = load_dynamics(dynamics_file)

    B = reshape(U_hat * B, 4, 256, 128)
    p1 = plot(1:256, 1:128, B[1,:,:]', title="Density Control")
    p2 = plot(1:256, 1:128, B[2,:,:]', title="X-Velocity Control")
    p3 = plot(1:256, 1:128, B[3,:,:]', title="Y-Velocity Control")
    p4 = plot(1:256, 1:128, B[4,:,:]', title="Energy Control")
    plot(p1,p2,p3,p4)
    savefig(output_img)
end

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
function continuous_prediction_error(A, B, U_hat, Ω, starting_points, T; verbose = true)
    n = size(U_hat, 1)
    q = size(Ω,1) - n
    average_err = Float64[]
    for s in starting_points
        verbose && println("Predicting from: ",s)
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

# Get the average prediction error over a run
average_continuous_prediction_error(A, B, U_hat, Ω, starting_points, T; verbose = true) = mean(continuous_prediction_error(A, B, U_hat, Ω, starting_points, T, verbose = verbose))

# Get the average continuous prediction err as a function of window size
function get_error_vs_window_size(train_dir, test_dir, train_window_sizes, test_starting_points, T; verbose = true)
    err = Float64[]
    modes = Int[]
    max_eig = Float64[]
    Ω, Xp = load_data(train_dir, 1:750)
    for sz in train_window_sizes
        verbose && println("training on window size: ", sz)
        A, B, phi, W, U_hat = DMDc(Ω[:,1:sz], Xp[:,1:sz], 0.99)
        push!(modes, size(A,1))
        push!(err, average_continuous_prediction_error(A,B,U_hat,Ω, test_starting_points, T, verbose = false))
        push!(max_eig, maximum(abs.(eigvals(A))))
    end
    err, modes, max_eig
end

# Get the average continuous prediction err, mode count and maximum eigenvalues as a function of window size
function plot_error_vs_window_size(train_dir, test_dir, train_window_sizes, test_starting_points, T, saveto)
    err, modes, max_eig = get_error_vs_window_size(train_dir, test_dir, train_window_sizes, test_starting_points, T)
    # Plot the results
    p1 = plot(train_window_sizes, err, label="", title="Average Prediction Error", xlabel="Window Size", ylabel="Average Prediction Error")

    p2 = plot(train_window_sizes, modes, label="", title="Number of Dynamic Modes vs. Window Size", xlabel="Window Size", ylabel = "Number of Dynamic Modes")

    p3 = plot(train_window_sizes, max_eig, label="", title="Maximum Eigenvalue of DM vs. Window Size", xlabel="Window Size", ylabel = "Maximum Eigenvalue")

    p = plot(p1, p2, p3)
    savefig(p, saveto)
end

# Plots the prediction accuracy over a specified range
function plot_prediction_accuracy(dynamics_file, dir, starting_points, T, output_img_name; read_control = true, data_index = Colon())
    # Load the dynamics
    print("Loading Dynamics...")
    A,B,U_hat = load_dynamics(dynamics_file)

    # Load the comparison data
    print("Loading Comparison Data...")
    Ω = load_data(dir, 1:maximum(starting_points) + T, read_control = read_control, data_index = data_index, get_next_frame = false)

    plot_prediction_accuracy(A, B, U_hat, Ω, starting_points, T, output_img_name)
end

# Plots the prediction accuracy over a specified range
function plot_prediction_accuracy(A, B, U_hat, Ω, starting_points, T, output_img_name)
    # Compute the average running error
    print("Running preditions:")
    avg_err = continuous_prediction_error(A, B, U_hat, Ω, starting_points, T)

    # Plot the results
    print("Plotting...")
    p = plot(starting_points, log.(avg_err), xlabel = "Iteration", ylabel="log(Average Error)", title = string("Average Prediction Error"))
    savefig(output_img_name)
    println("done!")
end

function make_vid_from_solution(dir, iter_range, data_index, data_name, output_img; fps = 25)
    anim = @animate for iteration in iter_range
        dict = h5_to_dict(get_filename(dir,iteration))
        sol_data = dict["sol_data"][data_index,:,:]
        plot(1:256, 1:128, sol_data', title = string(data_name, " at Iteration: ", iteration))
    end

    gif(anim, output_img, fps=25)
end




