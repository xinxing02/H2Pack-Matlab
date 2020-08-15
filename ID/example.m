%%  Numerically low-rank matrix example
%   two well-separated point clusters
X0 = 4 * (1 - 2 *rand(2000,3));
Y0 = bsxfun(@plus, [12, 0, 0], 4 * (1-2*rand(4000, 3)));

%   kernel block defined by K(x,y) = 1/|x-y|.
dist = pdist2(X0, Y0);
A = 1./dist;    


%%  ID with given rank. 
k = 2000;
f = 1.02;
[U, J] = ID(A, 'rank', k, f);

%   approximation error 
error = norm(A - U * A(J, :), 'fro') / norm(A, 'fro');

%   absolute maximum of entries in U
entry = max(abs(U(:)));

%   print result
fprintf('Relative approx. error: %.3E\n', error);
fprintf('Maximum abs. entry in U: %.3f\n\n', entry);


%%  ID with given error threshold. 
tol = 1e-4;
f = 1.2;
[U, J] = ID(A, 'tol', tol, f);

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
f = 1.2;
[P, J, E] = ID(A, 'rank', k, f);

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
