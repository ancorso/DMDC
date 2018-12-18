using LinearAlgebra
using StatsBase

function sorted_eig(M)
    D,V = eigen(M)
    i = sortperm(D, lt = (x,y) -> real(x) < real(y), rev = true)
    D = D[i]
    V = V[:, i]
    return D, V
end

function downsample(x, samples)
    n = size(x, 1)
    step = convert(Int, floor(n / samples))
    return x[1:step:n, :]
end

# Take the dynamic mode decomposition
function DMD(X, Xp, thresh)
    # take the singular value decomposition
    U, S, V = svd(X)
    r = length(S[S .> thresh])

    # Truncate the matrices
    Ur, Sr, Vr = U[:,1:r], diagm(0 => S[1:r]), V[:, 1:r]

    # Computed linearized truncated dynamics
    sinv = inv(Sr)
    Atil = Ur' * Xp * Vr * sinv
    D, W = sorted_eig(Atil)
    phi = Xp * Vr * sinv * W

    return Atil, phi, D
end

function compressed_DMD(X, Xp, samples, thresh)
    Y = downsample(X, samples)
    Yp = downsample(Xp, samples)

    # take the singular value decomposition
    U, S, V = svd(Y)
    r = length(S[S .> thresh])
    r = 13

    # Truncate the matrices
    Ur, Sr, Vr = U[:,1:r], diagm(0 => S[1:r]), V[:, 1:r]

    # Computed linearized truncated dynamics
    sinv = inv(Sr)
    Atil = Ur' * Yp * Vr * sinv
    D, W = sorted_eig(Atil)
    phiy = Yp * Vr * sinv * W
    phix = Xp * Vr * sinv * W
    return Atil, phix, D
end


# Take the dynamic mode decomposition with control
function DMDc(Omega, Xp, thresh)
    n = size(Xp, 1)
    q = size(Omega, 1) - n

    # Compute the singular value decomposition of Omega (the input space)
    U, S, V = svd(Omega)
    r = length(S[S .> thresh]) # this is r_tilde

    # Truncate the matrices
    U_til, S_til, V_til = U[:,1:r], diagm(0 =>S[1:r]), V[:, 1:r]

    # Compute this efficient SVD of the output space Xp
    U, S, V = svd(Xp)
    r = length(S[S .> thresh]) # this is the threshod for the output space and should be less than rtil (seems to be the case for the testcase tried)
    U_hat = U[:,1:r] # note that U_hat' * U_hat \approx I


    U_1 = U_til[1:n, :]
    U_2 = U_til[n+1:n+q, :]
    sinv = inv(S_til)

    approxA = U_hat' * Xp * V_til * sinv * U_1' * U_hat
    approxB = U_hat' * Xp * V_til * sinv * U_2'

    D, W = eigen(approxA)
    phi = (Xp * V_til * sinv) * (U_1' * U_hat * W)

    return approxA, approxB, phi, W, U_hat
end

