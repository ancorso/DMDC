using LinearAlgebra
using StatsBase

# Get the normalized singular values from the matrix
function normalized_singular_values(X)
    _, Σ, _ = svd(X)
    cumsum(Σ[2:end].^2) / sum(Σ[2:end].^2)
end

# Get the number of modes to keep for a specified amount of the energy to be maintained.
# Rounds up to the closest even number
function num_modes(Σ, retained_energy)
    r = findfirst(cumsum(Σ[2:end].^2) / sum(Σ[2:end].^2) .> retained_energy)
    Int(round((r+1)/2, RoundUp)*2)
end

# Compute the dynamic mode decomposition with control
function DMDc(Ω, Xp, retained_energy = 0.99)
    n = size(Xp, 1)
    q = size(Ω, 1) - n

    # Compute the singular value decomposition of Omega (the input space)
    U, Σ, V = svd(Ω)
    r = num_modes(Σ, retained_energy)

    # Truncate the matrices
    U_til, Σ_til, V_til = U[:,1:r], diagm(0 =>Σ[1:r]), V[:, 1:r]

    # Compute this efficient SVD of the output space Xp
    U, Σ, V = svd(Xp)
    rp = num_modes(Σ, retained_energy)
    U_hat = U[:,1:rp] # note that U_hat' * U_hat \approx I


    U_1 = U_til[1:n, :]
    U_2 = U_til[n+1:n+q, :]
    Σ_inv = inv(Σ_til)

    A = U_hat' * Xp * V_til * Σ_inv * U_1' * U_hat
    B = U_hat' * Xp * V_til * Σ_inv * U_2'

    D, W = eigen(A)
    ϕ = (Xp * V_til * Σ_inv) * (U_1' * U_hat * W)

    return A, B, ϕ, W, U_hat
end

# Downsample the rows of matrix x, to get "nsamps" evenly spaces samples
function downsample(x, nsamps)
    return x[range(1, size(x, 1), length=nsamps), :]
end

function compressed_DMD(X, Xp, samples, thresh)
    Y = downsample(X, samples)
    Yp = downsample(Xp, samples)

    # take the singular value decomposition
    U, Σ, V = svd(Y)
    r = length(Σ[Σ .> thresh])

    # Truncate the matrices
    Ur, Σr, Vr = U[:,1:r], diagm(0 => Σ[1:r]), V[:, 1:r]

    # Computed linearized truncated dynamics
    Σ_inv = inv(Σr)
    A = Ur' * Yp * Vr * Σ_inv
    D, W = eigen(A)
    ϕy = Yp * Vr * Σ_inv * W
    ϕx = Xp * Vr * Σ_inv * W
    return A, ϕx, D
end

