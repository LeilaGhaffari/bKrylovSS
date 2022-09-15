# Arnoldi process
function arnoldi(A, u, m; tol=1e-12)
    n = size(A, 1)
    Q = zeros(n, min(n, m + 1))
    H = zeros(min(n, m + 1), m)
    Q[:, 1] = u / norm(u)
    for j in 1:m
        v = A * Q[:, j]
        for i in 1:j
            H[i, j] = Q[:, i]' * v
            v -= Q[:, i] * H[i, j]
        end
        if norm(v) < tol
            break
        end
        H[j+1, j] = norm(v)
        Q[:, j+1] = v / H[j+1, j]
    end
    Q, H
end
