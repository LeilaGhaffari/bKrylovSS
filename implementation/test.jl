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
bArnoldi(A, B, M, L)
