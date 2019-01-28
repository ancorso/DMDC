include("load_data.jl")
include("dmdc.jl")

# Set window ranges
dir = "../cylinder_data/prop_control/"
data_index = Colon()
file_range = 1:100
thresh = 0.99 # Percent of energy retained in SVD

# Load in the data from the files
Ω, Xp = load_data(dir, file_range, true, data_index)

# Get the cumulative sum of energy from the singular values
s = normalized_singular_values(Ω)
sp = normalized_singular_values(Xp)
plot(s, title="Cumulative Singular Value (SV) Energy", ylabel="Cumulative Energy in SVs", xlabel="Number of SVs", label="Ω", xscale=:log10, markershape = :circ)
plot!(sp, label="X'", xscale=:log10, markershape = :circ)
savefig("dmdc_example_singular_values")

# perform the dynamic mode decomposition with control
A, B, phi, W, U_hat = DMDc(Ω, Xp, thresh)

# Consider the first four modes
r = size(A, 1)
dofs = data_index == Colon() ? 4 : 1
plots = []
for dof = 1:dofs
    for m = 1:r
        mode = reshape(real(phi[:,m]), 4, 256, 128)
        push!(plots, plot(1:256, 1:128, mode[dof,:,:]', title = string("Mode: ", m, " dof: ", dof)))
    end
end

# Plot and save the modes for each degree of freedom
plot(plots..., layout = (r, dofs), size = (600*dofs,400*r))
savefig("dmdc_example_modes")

