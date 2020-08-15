%%  Numerically low-rank matrix example
%   two well-separated point clusters
X0 = 4 * (1 - 2 *rand(1000,3));
Y0 = bsxfun(@plus, [12, 0, 0], 4 * (1-2*rand(1000, 3)));

%   kernel block defined by K(x,y) = 1/|x-y|.
a = 1;
eta = 1;
dim = 3;
kernel = @(coord)rpy_mex(coord, a, eta);
A = kernel({X0, Y0});

%%  IDdim with given rank. 
k = 100;
[U, J] = IDdim(A, dim, 'rank', k);

%   approximation error 
error = norm(A - U * A(J, :), 'fro') / norm(A, 'fro');

%   absolute maximum of entries in U
entry = max(abs(U(:)));

%   print result
fprintf('Relative approx. error: %.3E\n', error);
fprintf('Maximum abs. entry in U: %.3f\n\n', entry);

%%  ID with given error threshold. 
tol = 1e-4;
[U, J] = IDdim(A, dim, 'tol', tol);

%   approximation error 
error = norm(A - U * A(J, :), 'fro') / norm(A, 'fro');

%   obtained rank
k = length(J);

%   absolute maximum of entries in U
entry = max(abs(U(:)));

%   print result
fprintf('Approximation rank: %d\n', k);
fprintf('Relative approx. error: %.3E\n', error);
fprintf('Maximum abs entry in U: %.3f\n\n', entry);

%%  ID with given rank (with a different output)
k = 100;
[P, J, E] = IDdim(A, dim, 'rank', k);

%   construct U based on P and E (storing P and E is cheaper than storing U)
tmp = [eye(size(E, 2)); E];
U = tmp(P, :);

%   approximation error 
error = norm(A - U * A(J, :), 'fro') / norm(A, 'fro');

%   absolute maximum of entries in U
entry = max(abs(U(:)));

%   print result
fprintf('Relative approx. error: %.3E\n', error);
fprintf('Maximum abs. entry in U: %.3f\n\n', entry);
