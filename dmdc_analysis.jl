include("load_data.jl")
include("dmdc.jl")
using Statistics
using LinearAlgebra
using Plots; pyplot()
using ImageFiltering

# Load the dynamics file
function load_dynamics(file)
    f = h5open(file, "r")
    A, B, transform = read(f, "A"), read(f, "B"), read(f, "transform")
    close(f)

    A, B, transform
end

function plot_modes(dynamics_file, output_img, num_modes = nothing)
    _,_,transform = load_dynamics(dynamics_file)
    modes = pinv(transform)
    r = (num_modes == nothing) ? size(modes, 2) : num_modes
    dofs = 4
    plots = []
    for dof = 1:dofs
        for m = 1:r
            mode = reshape(real(modes[:,m]), 4, 256, 128)
            push!(plots, plot(1:256, 1:128, mode[dof,:,:]', title = string("Mode: ", m, " dof: ", dof), size=(600,400)))
        end
    end
    plot(plots..., layout = (dofs, r), size = (600*r,400*dofs))
    savefig(output_img)
end

function plot_B(dynamics_file, output_img)
    _,B,transform = load_dynamics(dynamics_file)

    B = reshape(pinv(transform) * B, 4, 256, 128)
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

function plot_control(dir, iterations, output_img_name)
    # Load in cost and control input at desired iterations
    control = get_scalar_sequences(dir, iterations, ["control_input"])

    # Plot the performance
    plot(iterations, control, title = "Control function vs. Iteration", xlabel="Iteration", ylabel = "Control")
    savefig(output_img_name)
end

function prediction_error(A, B, transform, Ω, s, T)
    n = size(transform, 2)
    q = size(Ω,1) - n
    x0 = Ω[1:n,s] # starting situation
    u = Ω[end,:]
    b = transform*x0
    detransform = pinv(transform)

    errl = Float64[]
    for i=s+1:s+T-1
        b = A*b + B*u[i+1]
        x = detransform * b

        push!(errl, norm(Ω[1:n, i] - x))
    end
    return errl
end


# A, B and transform are the dmdc params
# Ω is the true data
# starting_points are the iterations that the prediction should start at
# T is the prediction window
function continuous_prediction_error(A, B, transform, Ω, starting_points, T; verbose = true)
    println("size of A: ", size(A))
    average_err = Float64[]
    for s in starting_points
        verbose && println("Predicting from: ",s)
        push!(average_err, mean(prediction_error(A, B, transform, Ω, s, T)))
    end
    average_err
end

# Get the average prediction error over a run
average_continuous_prediction_error(A, B, transform, Ω, starting_points, T; verbose = true) = mean(continuous_prediction_error(A, B, transform, Ω, starting_points, T, verbose = verbose))

# Get the average continuous prediction err as a function of window size
function get_error_vs_window_size(train_dir, test_dir, train_window_sizes, test_starting_points, T; verbose = true)
    err = Float64[]
    modes = Int[]
    max_eig = Float64[]
    Ω, Xp = load_data(train_dir, 1:750)
    for sz in train_window_sizes
        verbose && println("training on window size: ", sz)
        A, B, phi, W, transform = DMDc(Ω[:,1:sz], Xp[:,1:sz], 0.99)
        push!(modes, size(A,1))
        push!(err, average_continuous_prediction_error(A, B, transform, Ω, test_starting_points, T, verbose = false))
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
    A, B, transform = load_dynamics(dynamics_file)

    # Load the comparison data
    print("Loading Comparison Data...")
    Ω = load_data(dir, 1:maximum(starting_points) + T, read_control = read_control, data_index = data_index, get_next_frame = false)

    plot_prediction_accuracy(A, B, transform, Ω, starting_points, T, output_img_name)
end

# Plots the prediction accuracy over a specified range
function plot_prediction_accuracy(A, B, transform, Ω, starting_points, T, output_img_name)
    # Compute the average running error
    print("Running predictions...")
    # avg_err_B = continuous_prediction_error(A, B, transform, Ω, starting_points, T)
    # avg_err_noB = continuous_prediction_error(A, zeros(size(B)), transform, Ω, starting_points, T)

    avg_err_B = prediction_error(A, B, transform, Ω, starting_points, T)
    avg_err_noB = prediction_error(A, zeros(size(B)), transform, Ω, starting_points, T)

    # Plot the results
    print("Plotting...")
    p = plot(log.(avg_err_B), xlabel = "Iteration", ylabel="log(Average Error)", title = string("Prediction Error from Initial State"), label = "With Control", legend = :bottomright)
    plot!(log.(avg_err_noB), label = "No Control")
    savefig(output_img_name)
    println("done!")
end


# Plots the continuous average prediction accuracy over a specified range
function plot_continuous_prediction_accuracy(dynamics_file, dir, starting_points, T, output_img_name; read_control = true, data_index = Colon())
    # Load the dynamics
    print("Loading Dynamics...")
    A, B, transform = load_dynamics(dynamics_file)

    # Load the comparison data
    print("Loading Comparison Data...")
    Ω = load_data(dir, 1:maximum(starting_points) + T, read_control = read_control, data_index = data_index, get_next_frame = false)

    plot_continuous_prediction_accuracy(A, B, transform, Ω, starting_points, T, output_img_name)
end

# Plots the moving average prediction accuracy
function plot_continuous_prediction_accuracy(A, B, transform, Ω, starting_points, T, output_img_name)
    # Compute the average running error
    print("Running predictions...")
    avg_err_B = continuous_prediction_error(A, B, transform, Ω, starting_points, T)
    avg_err_noB = continuous_prediction_error(A, zeros(size(B)), transform, Ω, starting_points, T)

    # Plot the results
    print("Plotting...")
    p = plot(log.(avg_err_B), xlabel = "Iteration", ylabel="log(Average Error)", title = string("Average Prediction Error over Next $T Timesteps"), label = "With Control", legend = :bottomright)
    plot!(log.(avg_err_noB), label = "No Control")
    savefig(output_img_name)
    println("done!")
end


function plot_h5(file, data_index, output_img, title)
    dict = h5_to_dict(file)
    sol_data = dict["sol_data"][data_index,:,:]
    sol_data = imfilter(sol_data, Kernel.gaussian(2))
    plot(1:256, 1:128, sol_data', title = title, xlabel="X", ylabel="Y")
    savefig(output_img)
end


function make_vid_from_solution(dir, iter_range, data_index, data_name, output_img; fps = 25)
    anim = @animate for iteration in iter_range
        dict = h5_to_dict(get_filename(dir,iteration))
        sol_data = dict["sol_data"][data_index,:,:]
        sol_data = imfilter(sol_data, Kernel.gaussian(2))
        plot(1:256, 1:128, sol_data', title = string(data_name, " at Iteration: ", iteration), xlabel="X", ylabel="Y")
    end

    gif(anim, output_img, fps=25)
end

function find_max_sol_file_num(dir, base, ext)
    files = readdir(dir)
    num = maximum([parse(Int, match(Regex("$(base)([0-9]{4})$(ext)"), f).captures[1]) for f in files])
end

function find_dynamics_file(dir)
    files = readdir(dir)
    for f in files
        m = match(r".*dynamics.*\.h5", f)
        if m != nothing
            return m.match
        end
    end
    return "NO MATCH FOUND"
end

if length(ARGS) > 0
    dir = "sol_data/"
    max_file_num = find_max_sol_file_num(dir, "sol_data_", ".h5")
    iter_range = 1:max_file_num
    T = 32
    dynamics_file = find_dynamics_file(".")
    println("Dynamics file: ", dynamics_file)


    for option in ARGS
        println("option processing: ", option)
        if option == "make_vid_from_solution"
            make_vid_from_solution(dir, iter_range, 3, "Y-Vel", "solution_vid.gif")
        elseif option == "plot_prediction_accuracy"
            plot_prediction_accuracy(dynamics_file, dir, 50, min(200,max_file_num), "prediction_accuracy")
        elseif option == "plot_continuous_prediction_accuracy"
            plot_continuous_prediction_accuracy(dynamics_file, dir, 50:3:max_file_num, T, "continuous_prediction_accuracy")
        elseif option == "plot_suppression_performance"
            plot_suppression_performance(dir, iter_range, "suppression_performance")
        elseif option == "plot_B"
            plot_B(dynamics_file, "control_response")
        else
            @error string("Unrecognized command: ", option)
        end
    end

end





