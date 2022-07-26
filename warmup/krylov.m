function km = krylov(n, M)

    lambda = 10 + (1:n);
    A = diag(lambda) + triu(rand(n), 1);
    b = rand(n, 1);
    
    km = b;
    for m = 1:M
        v = A * km(:, m);
        km(:, m+1) = v / norm(v);
    end    
end