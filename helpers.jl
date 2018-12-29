# find the point of max difference between timesteps and get the history at that point
function exact_sol_at_point(X)
    mx, index = findmax(X[:,1] - X[:,2])
    return index, X[index, :]
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

# Convert a global index into an ijk value for the given dimensions
function to_ijk(g_ind::Int, n_pts::Array{Int,1})
    dims = length(n_pts)
    inds = ones(Int, dims)
    g_ind -= 1
    slabs = [1, cumprod(n_pts)[1:end-1]...]
    for n = dims:-1:1
        inds[n] = floor(g_ind / slabs[n]) + 1
        g_ind = g_ind % slabs[n]
    end
    inds
end

# Convert an ijk value to an index
to_index(ijk::Array{Int,1}, n_pts::Array{Int,1}) = dot(ijk.-1, [1, cumprod(n_pts)[1:end-1]...]) .+ 1

# cap complex or real values to the specified magnitude
function cap(n, cap_val = 1)
    absn = abs(n)
    if absn > cap_val
        n /= absn / cap_val
    end
    return n
end