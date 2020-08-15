function Y = sample_densify(Y0, rho, domain)

N = size(Y0, 1);
dim = size(Y0, 2);

D = squareform(pdist(Y0)) + 1e10 * eye(N);
mindist = min(D, [], 2);

%   random selection around each point
tmp = 1 - 2*rand(2*rho*N, dim);
% %   normalize the selection (on a sphere)
% tmp = bsxfun(@rdivide, tmp, sum(tmp.^2, 2).^(1/2));
%   set up the distance
tmp = bsxfun(@times, tmp, kron(mindist, 0.33*ones(2*rho, 1)));
%   move to the corrsponding center
tmp = bsxfun(@plus, tmp, kron(Y0, ones(2*rho, 1)));
%   eliminate the points outside the domain
tmp(domain(tmp) <= 0, :) = [];

Y = [Y0; tmp(randperm(size(tmp,1), min(rho*N, size(tmp,1))), :)];

end