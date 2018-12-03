using HDF5

# Construct the h5 filename at the desired iteration for the given naming and storage scheme
function get_filename(iteration, dir, base = "sol_data_", pad = 4, ext = ".h5")
    string(dir, base, lpad(string(iteration), 4, '0'), ext)
end

# Construct the h5 filename of the controlled cylinder at the desired iteration
get_controlled_cyl_filename(iteration) = get_filename(iteration, "data/rot_cyl_re50/")
# Construct the h5 filename of the static cylinder at the desired iteration
get_static_cyl_filename(iteration) = get_filename(iteration, "data/static_cyl/")
# Construct the h5 filename of the dfc run
get_dfc_test_filename(iteration) = get_filename(iteration, "dfc_run_data/")

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
function fill_control_snapshots(Omega, Xp)
    l, n = size(Omega, 2), size(Xp, 1)
    for i=1:l+1
        u,x = read_snapshot(get_controlled_cyl_filename(i))
        (i <= l) && (Omega[1:n,i] = x[:])
        (i <= l) && (Omega[n+1:end,i] = [u])
        (i > 2) && (Xp[:,i-1] = x[:])
    end
end

# Construct the solution data matrices for the static case
function fill_static_snapshots(X, Xp)
    n, last = size(X)
    for i=0:last
        x = read_snapshot(get_static_cyl_filename(i))[1]
        (i < last) && (X[:,i+1] = x[:])
        (i > 0) && (Xp[:,i] = x[:])
    end
end



