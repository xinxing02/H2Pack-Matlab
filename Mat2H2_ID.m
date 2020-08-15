function h2mat = Mat2H2_ID(A, htree, type, par, alpha, kdim)

%   Construction of the symmetric H^2 approximation with ID. 
%
%   Hierarchical Structure: htree
%
%   Compression Paremeter for ID: type, par;
%
%   Admissible condition: alpha
%
%   Dimension of the kernel : kdim
%
%   Output: ID for eahc HSS block row: A(, :) = P [eye; E] A(I, :)
%           near, k * 2 array, list all the dense blocks. 
%                   (i,j) denotes cluster i & cluster j
%           far,  k * 2 array, list all the admissible blocks. 
%

%   check input
if nargin < 6
    kdim = 1;
end

%   Basic info
parent   = htree.parent;
children = htree.children;
level    = htree.level;
nodelvl  = htree.nodelvl;
leafnode = htree.leafnode;
nlevel   = length(level);
if isfield(htree, 'mcluster')
    mcluster = htree.mcluster;
else
    mcluster = [kdim * htree.cluster(:,1) - (kdim-1), kdim * htree.cluster(:,2)];
end

%   Initialization
tnode = length(parent);
D = cell(tnode);
B = cell(tnode);
U = cell(tnode,1); 
I = cell(tnode,1); 
if (isscalar(par))
    par = par * ones(nlevel,1);
end


%   Inadmissible and admissible blocks
[near, far] = H2__block_partition_box(htree, alpha);
minlvl = min(nodelvl(far(:)));

%   Admissible node list at each level's construction.
admnode = H2__admissible_nodelist(htree, near, far);


%   Initialize the row indices for leafnode clusters
for i = 1 : length(leafnode)
    node = leafnode(i);
    I{node} = mcluster(node, 1) : mcluster(node, 2);
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
        
        %   Column indices for the compression
        if isempty(admnode{node})
            continue;
        else
            colidx = horzcat(I{admnode{node}});            
        end    
        
        %   Randomized method 
        tmpA = A(I{node}, colidx) * randn(length(colidx), length(I{node})); 
        colnorm = sqrt(sum(tmpA.^2, 1));
        tmpA = bsxfun(@times, tmpA, 1./colnorm);        

        %   Row ID for the block rows.
        if kdim == 1
            [U{node}, subidx] = ID(tmpA, type, par(i));
        else
            [U{node}, subidx] = IDdim(tmpA, kdim, type, par(i));            
        end
        I{node} = I{node}(subidx);
    end
end

%%   Inadmissible/Intermediate blocks
%   diagonal blocks   
for i = 1 : length(leafnode)
    node = leafnode(i);
    idx  = mcluster(node,1) : mcluster(node,2);
    D{node, node} = A(idx, idx);
end

%   off-diagonal in-admissible blocks
for i = 1 : size(near, 1)
    c1 = near(i,1);
    c2 = near(i,2);
    idx1 = mcluster(c1, 1): mcluster(c1, 2);
    idx2 = mcluster(c2, 1): mcluster(c2, 2);
    D{c1, c2} = A(idx1, idx2);
end   

%   intermediate blocks
for i = 1 : size(far, 1)
    c1 = far(i,1);
    c2 = far(i,2);
    idx1 = mcluster(c1, 1): mcluster(c1, 2);
    idx2 = mcluster(c2, 1): mcluster(c2, 2);    
    if nodelvl(c1) == nodelvl(c2)
        %   compression on both sides
        B{c1,c2} = A(I{c1}, I{c2});
    elseif nodelvl(c1) > nodelvl(c2)  
        %   c2 is the leafnode at higher level. only compress on c1's side
        B{c1,c2} = A(I{c1}, idx2);
    else
        %   c1 is the leafnode at higher level. only compress on c2's side
        B{c1,c2} = A(idx1, I{c2});
    end
end

%%  Wrapup
h2mat.D = D;
h2mat.B = B;
h2mat.U = U;
h2mat.I = [];
h2mat.near = near;
h2mat.far = far;
h2mat.minlvl = minlvl;
h2mat.JIT = false;
h2mat.alpha = alpha;
h2mat.type = type;
h2mat.par = par;
h2mat.kernel = nan;
h2mat.kdim = kdim;
h2mat.storage = H2__storage_cost(h2mat, htree);
h2mat.rankinfo = H2__rank(h2mat, htree);

end