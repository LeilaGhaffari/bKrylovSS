# Reference:
#   Iterative Methods for Sparse Linear Systems
#     Yousef Saad
#     page 219 -> TODO: implement BlockArnoldi–Ruhe’s variant (ALGORITHM 6.24, page 220)

# Helper function for updating the Arnoldi matrices
iter(i, p) = (i-1)*p+1:i*p

# Block Arnoldi
function bArnoldi(A, X0, m; tol=1e-12)
    n = size(A, 1)
    p = size(X0, 2)
    W = zeros(n, p * (m + 1))   # W = [V1, ..., Vm+1]
    H = zeros(p * (m + 1), p * m) # H = [S1, ..., Sm]
    f = qr(X0)
    W[:, 1:p] = f.Q[:, 1:p] # GMRES will need to keep f.R
    S0 = f.R
    for j in 1:m
        iterJ = iter(j, p)
        zz = A * W[:, iterJ]
        for i in 1:j
            iterI = iter(i, p)
            H[iterI, iterJ] = W[:, iterI]' * zz
            zz -= W[:, iterI] * H[iterI, iterJ]
        end
        iterJp1 = iter(j + 1, p)
        V, S = qr(zz)
        if norm(S) < tol
            break
        end
        W[:, iterJp1] = V[:, 1:p]
        H[iterJp1, iterJ] = S
    end
    W, H, S0
end
