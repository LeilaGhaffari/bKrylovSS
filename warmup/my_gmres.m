clc
clear all
maxit = 120;
rtol = 1e-8;
d = 50;
A = d^2 * gallery('poisson', d);
n = size(A, 1)
b = ones(n, 1);
rest = [maxit 20 40 60];
for j = 1:4
    [~, ~, ~, ~, rv] = gmres(A, b, rest(j), rtol, maxit/rest(j));
    semilogy(0:length(rv)-1, rv, '-'), hold on
end