function h2mat = Mat2H2_ID_MatrixFree(kernel, htree, type, par, alpha, JIT_flag)

%   Construction of the symmetric H^2 approximation with ID. 
%
%   Hierarchical Structure: htree
%
%   Compression Paremeter for ID: type, par;
%
%   Admissible condition: alpha
%
%   Output: ID for eahc HSS block row: A(, :) = P [eye; E] A(I, :)
%           near, k * 2 array, list all the dense blocks. 
%                   (i,j) denotes cluster i & cluster j
%           far,  k * 2 array, list all the admissible blocks. 

if nargin < 6
    JIT_flag =false;
end

%   Basic info
parent   = htree.parent;
children = htree.children;
level    = htree.level;
nodelvl  = htree.nodelvl;
leafnode = htree.leafnode;
nlevel   = length(level);
cluster  = htree.cluster;
enbox    = htree.enbox;
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

%   Initialize the row indices for leafnode clusters
for i = 1 : length(leafnode)
    node = leafnode(i);
    I{node} = cluster(node, 1) : cluster(node, 2);
end

%   Inadmissible and admissible blocks
[near, far] = H2__block_partition_box(htree, alpha);
minlvl = min(nodelvl(far(:)));

%   Admissible node list at each level's construction.
admnode = H2__admissible_nodelist(htree, near, far);


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
        
        %   FF nodes indices for the compression
        if isempty(admnode{node})
            continue;
        else
            colidx = horzcat(I{admnode{node}});            
        end    
        Aobj = kernel({coord(I{node}, :), coord(colidx, :)});

        %   Randomized method 
        tmpA = Aobj * randn( size(Aobj,2), size(Aobj,1) ); 
        colnorm = sqrt(sum(tmpA.^2, 1));
        tmpA = bsxfun(@times, tmpA, 1./colnorm);  
        
        %   Row ID for the block rows.
        if kdim == 1
            [U{node}, subidx] = ID(tmpA, type, par(i));
        else
            [U{node}, subidx] = IDdim(tmpA, kdim, type, par(i));
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
h2mat.D = D;
h2mat.B = B;
h2mat.U = U;
h2mat.I = I;
h2mat.near = near;
h2mat.far = far;
h2mat.minlvl = minlvl;
h2mat.JIT = JIT_flag;
h2mat.alpha = alpha;
h2mat.type = type;
h2mat.par = par;
h2mat.kernel = kernel;
h2mat.kdim = kdim;
h2mat.storage = H2__storage_cost(h2mat, htree);
h2mat.rankinfo = H2__rank(h2mat, htree);
end