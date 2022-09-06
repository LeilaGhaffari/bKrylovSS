using LinearAlgebra
using SparseArrays

include("arnoldi.jl")
include("bArnoldi.jl")


# Set up the Shifted Linear Systems
n = 6
M = 3
L = 2
#^^^^^^^^^^^^^^^^^^^^^^^#
#  A*X_σ - X_σ * D = B  #
#_______________________#
# LHS
# -- A
λ = [i+10 for i in 1:n]
a = randn(n, n)
A = UpperTriangular(a) - diagm(diag(a)) + diagm(λ)
# -- Shift
σ = rand(1:9, L)
D = diagm(σ) 
# RHS
B = randn(n, L)

# Initialization
X_σ     = zeros(n, L) 
X_exact = zeros(n, L)
Res     = zeros(M, L)

# Test block Arnoldi
W, H = bArnoldi(A, B, M, L)

for m in 1:M
    s = zeros(L*(m+1), L)
    for i=1:L
        s[1, i] = norm(B[:, i])
    end
    z = H[1:L*(m+1), 1:L*m] \ s
    x = W[:, 1:L*m] * z
    for i=1:L
        Res[m, i] = norm(B[:, i] - A*x[:, i])
    end
    X_σ[:, :] = x
end

@show X_σ
