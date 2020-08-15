function Yp = H2__Proxy_Point_Surface(htree, dim, alpha, N)
%
%   Construct the proxy surface points for boxes at each level of htree.  
%   
%   alpha decides where the proxy surface lies
%
    
    %   correction for negative alpha which is for HSS
    alpha = max(1e-3, alpha);
     
    %   basic info
    level    = htree.level;
    enbox    = htree.enbox;
    nlevel   = length(level);
    
    %   uniform proxy points on the surface of the cubic box [-1, 1]^dim 
    Ysurf = gridpoint_on_boxsurface(-1*ones(1, dim), eye(dim), 2*ones(1,dim), N, dim);

    %   properly scale Ysurf for different levels
    Yp = cell(nlevel, 1);
    for i = 2 : nlevel
        %   half edge length for box at level i
        L1 = enbox{level{i}(1)}(2,1) / 2;
        %   half edge length of the proxy surface at level i
        L2 = (1 + 2*alpha) * L1;
        %   scale
        Yp{i} = Ysurf * L2;
    end
end
