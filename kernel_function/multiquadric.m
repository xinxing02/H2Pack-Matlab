function M = multiquadric(coord, alpha, scale)
%
%   Quadric Function F(x,y) = (1 + scale * (x-y)^2) ^ (-alpha)
%

if (nargin < 3)
    scale = 1;
end

if isnumeric(coord)   
    N = size(coord, 1);
    if N == 1
        M = 1;
        return ;
    end
    dist = pdist(coord);
    M = squareform((1 + scale * dist.^2).^(-alpha)) + eye(N);
elseif iscell(coord)&&length(coord)==2      
    M = pdist2(coord{1}, coord{2});
    M = (1 + scale * M.^2).^(-alpha);
end