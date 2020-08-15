function M = matern(coord, e)
%
%   Exponential Kernel matrix, K(x,y) = ?1+e*dist) * exp(-eps * ||x-y||^2)
%
if isnumeric(coord)
    N = size(coord, 1);
    if N == 1
        M = 1;
        return ;
    end
    dist = pdist(coord);
    M = squareform( (1+e*dist) .* exp(- e * dist) ) + eye(N);
elseif iscell(coord)&&length(coord)==2      
    M = pdist2(coord{1}, coord{2});
    M = (1 + e*M) .* exp(- e * M);
end
end