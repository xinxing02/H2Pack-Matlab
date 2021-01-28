function M = logfun(coord)
%
%   Exponential Kernel matrix, K(x,y) = exp(-eps * ||x-y||^2)
%
if isnumeric(coord)
    N = size(coord, 1);
    if N == 1
        M = 0;
        return ;
    end 
    dist = pdist(coord);
    M = squareform(log(dist));
    M(abs(M) > 1e20) = 0;
elseif iscell(coord)&&length(coord)==2      
    M = pdist2(coord{1}, coord{2});
    M = log(M);
    M(abs(M) > 1e20) = 0;
end
end