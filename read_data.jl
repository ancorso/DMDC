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

# Construct the solution data and control matrices
function fill_control_snapshots(Omega, Xp, dir)
    l, n = size(Omega, 2), size(Xp, 1)
    for i=1:l+1
        u,x = read_snapshot(get_filename(i, dir))
        (i <= l) && (Omega[1:n,i] = x[:])
        (i <= l) && (Omega[n+1:end,i] = [u])
        (i > 2) && (Xp[:,i-1] = x[:])
    end
end

# Construct the solution data matrices for the static case
function fill_static_snapshots(X, Xp, dir)
    n, last = size(X)
    for i=0:last
        x = read_snapshot(get_filename(i, dir))[1]
        (i < last) && (X[:,i+1] = x[:])
        (i > 0) && (Xp[:,i] = x[:])
    end
end



