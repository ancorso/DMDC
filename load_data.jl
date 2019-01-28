using HDF5

# Construct the h5 filename at the desired iteration for the given naming and storage scheme
function get_filename(dir, iteration, base = "sol_data_", pad = 4, ext = ".h5")
    string(dir, base, lpad(string(iteration), 4, '0'), ext)
end

# Read an h5 file and conver it to a dictionary
function h5_to_dict(filename)
    f = h5open(filename, "r")
    output = Dict()
    for n in names(f)
        output[n] = read(f, n)
    end
    return output
end

# Get the control sequence from a set of data files store in dir
function get_scalar_sequences(dir, iterations, names::Array{String,1})
    vecs = [Array{Float64}(undef, length(iterations)) for n in names]
    i = 1
    for file_iter in iterations
        f = h5open(get_filename(dir, file_iter), "r")
        for n in length(names)
            vec[n][i] = read(f, names[n])
        end
        close(f)
        i = i+1
    end
    vecs
end

# Load data from the circular cylinder with or without control
function load_data(dir, iterations; data_index = Colon(), get_next_frame = true, read_control = true)
    next_iterations = iterations .+ 1 # Get the indices of the next frame
    nc = length(iterations) # The number of columns in the output matrices
    d = h5_to_dict(get_filename(dir, iterations[1])) # read the first data to get the size of the output
    n = length(d["sol_data"]) # The number of rows in Xp
    nr =  read_control ? n + length(d["control_input"]) : n # The number of rows in X
    X = Array{Float64,2}(undef, nr, nc)
    Xp = get_next_frame ?  Array{Float64,2}(undef, n, nc) : Float64[]

    for i in unique([iterations...,next_iterations...])
        d = h5_to_dict(get_filename(dir, i))
        xi = findall(iterations .== i)
        xpi = findall(next_iterations .== i)
        if !isempty(xi)
            X[1:n, xi[1]] = d["sol_data"][data_index, :, :][:]
            (read_control) && (X[n+1:nr, xi[1]] = [d["control_input"]])
        end
        (get_next_frame && !isempty(xpi)) && (Xp[:, xpi[1]] = d["sol_data"][data_index, :, :][:])
    end
    get_next_frame ? (X, Xp) : X
end


