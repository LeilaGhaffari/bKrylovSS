using LinearAlgebra
using SparseArrays
using Plots

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

# Helper function for updating the Arnoldi matrices
iter(i, p) = (i-1)*p+1:i*p

# Block Arnoldi based on 
#   Iterative Methods for Sparse Linear Systems, Yousef Saad, section 6.12
function bArnoldi(A, X0, m, p)
    n = size(A, 1)
    U = zeros(n, p*(m+1))
    H = zeros(p*(m+1), p*m)
    for i=1:p
        U[:, i] = X0[:, i] / norm(X0[:, i])
    end
    for j in 1:m
        iterJ = iter(j, p)
        W = A * U[:, iterJ]
        for i in 1:j
            iterI = iter(i, p)
            H[iterI, iterJ] = U[:, iterI]' * W
            W -= U[:, iterI] * H[iterI, iterJ]
        end
        iterJp1 = iter(j+1, p)
        Q, R = qr(W)
        U[:, iterJp1] = Q[:, 1:p]
        H[iterJp1, iterJ] = R
    end
    U, H
end

# Test block Arnoldi
bArnoldi(A, B, M, L)
