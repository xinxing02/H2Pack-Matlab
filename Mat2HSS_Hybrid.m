function hssmat = Mat2HSS_Hybrid(kernel, htree, Yp, type, par, alpha_pp, JIT_flag)

%   Construction of the symmetric H^2 approximation with ID. 
%
%   Kernel Info: kernel, coord
%
%   Hierarchical Structure: htree
%
%   Compression Paremeter for ID: type, par;
%
%   Admissible condition for far-field: alpha
%
%   Output: ID for eahc HSS block row: A(, :) = P [eye; E] A(I, :)
%           near, k * 2 array, list all the dense blocks. 
%                   (i,j) denotes cluster i & cluster j
%           far,  k * 2 array, list all the admissible blocks. 

if nargin < 7
    JIT_flag =false;
end

%   Basic info
parent   = htree.parent;
children = htree.children;
level    = htree.level;
nodelvl  = htree.nodelvl;
leafnode = htree.leafnode;
enbox    = htree.enbox;
cluster  = htree.cluster;
nlevel   = length(level);
coord    = htree.coord;

%   Kernel dimension
ptdim = size(enbox{1}, 2);
kdim = size(kernel(randn(2, ptdim)), 1) / 2;

%   Initialization
tnode = length(parent);
D = cell(tnode);
B = cell(tnode);
U = cell(tnode,1); 
I = cell(tnode,1); 
if (isscalar(par))
    par = par * ones(nlevel,1);
end


%   NN and FF blocks (used for the hybrid construction)
[near0, far0] = H2__block_partition_box(htree, alpha_pp);

%   NN and FF (used for the HSS matrix representation)
[near, far] = H2__block_partition_box(htree, 0);
minlvl = min(nodelvl(far(:)));

%   Inadmissible node list at each level's construction.
inadmnode = H2__inadmissible_nodelist(htree, near0, far0);

%   Initialize the row indices for leafnode clusters
for i = 1 : length(leafnode)
    node = leafnode(i);
    I{node} = cluster(node, 1) : cluster(node, 2);
end

%%   Hierachical construction level by level
for i = nlevel : -1 : minlvl
    
    %   Update row indices I{node} associated with the clusters at ith level
    for j = 1 : length(level{i})
        node = level{i}(j);
        child_node = children(node, ~isnan(children(node,:)));
        if ~isempty(child_node) 
            I{node} = horzcat(I{child_node});
        end 
    end
    
    %   Compression at the ith level
    for j = 1 : length(level{i})
        node = level{i}(j);
        boxcenter = enbox{node}(1,:) + 0.5 * enbox{node}(2,:);
        tmp_coord = bsxfun(@minus, coord(I{node},:), boxcenter);
        
        %   proxy point matrix
        if i > 2
            A_sam = kernel({tmp_coord, Yp{i}});   
        else
            A_sam = zeros(size(tmp_coord,1), 0);    % for second level, no need for proxy points
        end
        
        %   near field blocks
        if ~isempty(inadmnode{node})
            near_idx = horzcat(I{inadmnode{node}});
            Ynear = bsxfun(@minus, coord(near_idx, :), boxcenter);
            A_nf = kernel({tmp_coord, Ynear});
        else
            A_nf = zeros(kdim*size(tmp_coord, 1), 0);
        end
        
        %   ID compression
        A_hybrid = [A_sam, A_nf];
        
        %   Randomization
        A_hybrid = A_hybrid * randn(size(A_hybrid,2), size(A_hybrid,1));
        colnorm = sqrt(sum(A_hybrid.^2, 1));
        A_hybrid = bsxfun(@times, A_hybrid, 1./colnorm);        
        
        if kdim == 1
            [U{node}, subidx] = ID(A_hybrid, type, par(i));
        else
            [U{node}, subidx] = IDdim(A_hybrid, kdim, type, par(i));
            subidx = subidx(kdim:kdim:end)/kdim;
        end
        I{node} = I{node}(subidx);                                            
    end
end

%%   Inadmissible/Intermediate blocks
if JIT_flag == false
    %   diagonal blocks   
    for i = 1 : length(leafnode)
        node = leafnode(i);
        idx  = cluster(node,1) : cluster(node,2);
        D{node, node} = kernel(coord(idx,:));
    end

    %   off-diagonal in-admissible blocks
    for i = 1 : size(near, 1)
        c1 = near(i,1);
        c2 = near(i,2);
        idx1 = cluster(c1, 1): cluster(c1, 2);
        idx2 = cluster(c2, 1): cluster(c2, 2);
        D{c1, c2} = kernel({coord(idx1, :), coord(idx2, :)});
    end   

    %   intermediate blocks
    for i = 1 : size(far, 1)
        c1 = far(i,1);
        c2 = far(i,2);
        idx1 = cluster(c1, 1): cluster(c1, 2);
        idx2 = cluster(c2, 1): cluster(c2, 2);    
        if nodelvl(c1) == nodelvl(c2)
            %   compression on both sides
            B{c1,c2} = kernel({coord(I{c1},:), coord(I{c2},:)}); 
        elseif nodelvl(c1) > nodelvl(c2)  
            %   c2 is the leafnode at higher level. only compress on c1's side
            B{c1,c2} = kernel({coord(I{c1},:), coord(idx2, :)});
        else
            %   c1 is the leafnode at higher level. only compress on c2's side
            B{c1,c2} = kernel({coord(idx1, :), coord(I{c2}, :)});
        end
    end
end

%%  Wrapup
hssmat.D = D;
hssmat.B = B;
hssmat.U = U;
hssmat.I = I;
hssmat.near = near;
hssmat.far = far;
hssmat.minlvl = minlvl;
hssmat.JIT = JIT_flag;
hssmat.alpha = 0;
hssmat.alpha_pp = alpha_pp;
hssmat.type = type;
hssmat.par = par;
hssmat.kernel = kernel;
hssmat.kdim = kdim;
hssmat.Yp = Yp;
hssmat.storage = H2__storage_cost(hssmat, htree);
hssmat.rankinfo = H2__rank(hssmat, htree);
end
