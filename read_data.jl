using HDF5

# Construct the h5 filename at the desired iteration for the given naming and storage scheme
function get_filename(iteration, dir, base = "sol_data_", pad = 4, ext = ".h5")
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

# Read the h5 file at the specified iteration and return the solution data and the control input
function read_snapshot(filename)
    # Note that in each hdf5 file, there is density, ux, uy, energy
    f = h5open(filename, "r")
    output = []
    for n in names(f)
        push!(output, read(f, n))
    end
    return output
end

function get_control_sequence(iter_rng, dir)
    control = []
    for iteration in iter_rng
        println("Loading file: ", iteration)
        dict = h5_to_dict(get_filename(iteration, dir))
        push!(control, dict["control_input"])
    end
    control
end

# Construct the solution data and control matrices
function fill_control_snapshots(Omega, Xp, dir, data_index = Colon(), ts = 1, te = nothing)
    l, n = size(Omega, 2), size(Xp, 1)
    (te == nothing) && (te = l+1)
    for iter=ts:te
        data_dict = h5_to_dict(get_filename(iter, dir))
        i = iter - ts + 1
        u = data_dict["control_input"]
        x = data_dict["sol_data"]
        (i <= l) && (Omega[1:n,i] = x[data_index, :, :][:])
        (i <= l) && (Omega[n+1:end,i] = [u])
        (i > 2) && (Xp[:,i-1] = x[data_index, :, :][:])
    end
end

# Construct the solution data matrices for the static case
function fill_static_snapshots(X, Xp, dir, index = Colon())
    n, last = size(X)
    for i=0:last
        x = read_snapshot(get_filename(i, dir))[1]
        (i < last) && (X[:,i+1] = x[index, :, :][:])
        (i > 0) && (Xp[:,i] = x[index, :, :][:])
    end
end



