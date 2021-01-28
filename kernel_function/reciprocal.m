function M = reciprocal(coord, alpha)
%
%   Reciprocal Kernel matrix, K(x,y) = 1 / ||x-y||^alpha
%
if isnumeric(coord)
    N = size(coord, 1);
    if N == 1
        M = 0;
        return ;
    end
    dist = pdist(coord);
    M = squareform(dist.^(-alpha));
    M(M > 1e20) = 0;
elseif iscell(coord)&&length(coord)==2      
    M = pdist2(coord{1}, coord{2});    
    M = 1./(M.^alpha);
    M(M > 1e20) = 0;
end
end