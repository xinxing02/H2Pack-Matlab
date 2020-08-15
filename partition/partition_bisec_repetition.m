function [parent, children, cluster, enbox, permutation, n] = ...
                        partition_bisec_repetition(coord, minSize, dim, box)

%   Hierarchical partitioning of the given points 'coord' in dim-D space.
%   Post-order tree is used
%
%   Input: 
%       coord,      N*3 matrix
%       minSize,    minimum size for the leaf nodes
%       dim,        dimension of the points. 
%       box,        2*dim matrix representing the box that encloses the
%                   point set 'coord'. 
%                   box(1, :) gives the coordinates of the corner of the box 
%                   with minimal x/y/z. 
%                   box(2, :) gives the edge length of the box along each
%                   dimension. 
%   
%   Output:
%       n:          number of nodes in the tree
%       parent:     n*1 vector, indicating the parent of the node. 
%       children:   n*2^dim matrix, indicating the children of the node.
%       cluster:    n*2 matrix, defining the indices at each cluster,
%                       [start-point, end-point]
%       enbox:      n*1 cell, each contains the box data that encloses the
%                       the associated cluster
%       permutation:n*1 vector, the permutation of the coordinates that makes 
%                       the clusters in order 

%   Check dimension
if size(coord,2) < dim
    disp 'Inconsistent coordinate dimension with the given dim'
    parent = []; children = []; cluster = []; enbox = {}; permutation = []; n = 0;
    return ;
else
    coord(:,(dim+1):end) = [];
end

%   clustering the unique coordinate first
[unique_coord, ~, ic] = unique(coord, 'rows');
if nargin == 3
    [parent, children, cluster, enbox, permutation, n] = ...
                        partition_bisec(unique_coord, minSize, dim);
else
    [parent, children, cluster, enbox, permutation, n] = ...
                        partition_bisec(unique_coord, minSize, dim, box);
end
                    
%   modify the cluster & permutation
iperm = invertperm(permutation);
newperm = zeros(size(coord, 1), 1);
newcluster = zeros(size(cluster,1),2);
idx = 0;
for i = 1 : size(cluster, 1)
    %   i is a leaf node
    if isnan(children(i, 1))
        newcluster(i, 1) = idx + 1;
        for j = cluster(i,1) : cluster(i, 2)
            rep_idx = find(ic == iperm(j));
            newperm(rep_idx) = (idx+1):(idx+length(rep_idx));
            idx = idx + length(rep_idx);        
        end
        newcluster(i, 2) = idx;
    %   i is a non-leaf node
    else
        nchild = find(~isnan(children(i,:)), 1, 'last');
        newcluster(i, 1) = newcluster(children(i,1), 1);
        newcluster(i, 2) = newcluster(children(i,nchild), 2);
    end
end

permutation = newperm;
cluster     = newcluster;
end







