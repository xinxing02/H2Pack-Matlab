function [near, far] = H2__block_partition_box(htree, alpha)
%
%   Calculate the node pairs [I,J] s.t. box{I} and box{J} 
%   are admissible or inadmissible pair in the H2 matrix representation. 
%
%   Input:
%       htree:      information of hierarchical partitioning
%       alpha:      admissibility coeff. 
%
%   Output: 
%       near,   N * 2 array with each row [I,J] denotes one inadmissible 
%               node pair.
%       far,    N * 2 array with each row [I,J] denotes one admissible 
%               node pair.
%
%
%   For any (i,j) in the far list, it must be one of the following case
%       1. if i and j are at the same level, then parent of i and parent of
%       j are not admissible. 
%       2. if i and j are at different level, then at least one of them 
%       must be leafnode. Furthermore, the level of the leafnode is not
%       less than that of the other and the leafnode is not admissible to
%       the parent of the other. 
%       

enbox = htree.enbox;
parent = htree.parent;
children = htree.children;
root = length(parent);
[near, far] = H2SelfIntersection(root);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [near, far] = H2Intersection(p1, p2)
%   Find the list of admissible/inadmissible blocks for two sub-trees
%   rooted at node1/node2.

%   initialization
near = zeros(0, 2);
far = zeros(0, 2);

%   admissible pair
if (isadmissible_box(enbox{p1}, enbox{p2}, alpha))
    near = [];
    far = [p1, p2];
    return ;
end

%   two leaf node and not admissible
if (isnan(children(p1, 1)) && isnan(children(p2, 1)))
    near = [p1, p2];
    far = [];
    return ;
end

%   p1 leaf and p2 nonleaf
if (isnan(children(p1, 1)) && ~isnan(children(p2, 1)))
    nchild = find(~isnan(children(p2,:)), 1, 'last');
    for i = 1 : nchild
        [nn, ff] = H2Intersection(p1, children(p2, i));
        near = vertcat(near, nn);
        far  = vertcat(far, ff);
    end
    return ;
end

%   p1 nonleaf and p2 leaf
if (~isnan(children(p1, 1)) && isnan(children(p2,1)))
    nchild = find(~isnan(children(p1,:)), 1, 'last');
    for i = 1 : nchild
        [nn, ff] = H2Intersection(children(p1, i), p2);
        near = vertcat(near, nn);
        far  = vertcat(far, ff);
    end
    return ;
end

%   both p1 and p2 nonleaf 
if (~isnan(children(p1, 1)) && ~isnan(children(p2,1)))
    nchild1 = find(~isnan(children(p1,:)), 1, 'last');
    nchild2 = find(~isnan(children(p2,:)), 1, 'last');
    for i = 1 : nchild1
        for j = 1 : nchild2
            [nn, ff] = H2Intersection(children(p1, i), children(p2, j));
            near = vertcat(near, nn);
            far  = vertcat(far, ff);
        end
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [near, far] = H2SelfIntersection(p)
%   Find the list of admissible/inadmissible blocks within the sub-trees
%   rooted at p

%   initialization
near = zeros(0, 2);
far = zeros(0, 2);

%   leaf node
if isnan(children(p, 1))
    near = [];
    far = [];
    return;
end

%   nonleaf node
nchild = find(~isnan(children(p,:)), 1, 'last');

%   self-intersection
for i = 1 : nchild
    [nn, ff]  = H2SelfIntersection(children(p,i));
    near = vertcat(near, nn);
    far  = vertcat(far, ff);
end

%   inter-interaction
for i = 1 : nchild
    for j = (i+1) : nchild
        [nn, ff]  = H2Intersection(children(p,i), children(p,j));
        near = vertcat(near, nn);
        far  = vertcat(far, ff);
    end
end
end %End for H2SelfIntersection()

end %End for H2block()

