function [NFList, FFList] =   H2__interaction_list(htree, row_htree, col_htree, alpha)
%
%   Given a partition tree and two sets of leaf nodes, row_leafnode, and
%   col_leafnode, construct the interaction list for the matrix block
%   A(row_leafnode, col_leafnode), or equivalently the block list. 
%
%   note: no symmetry is available. For any constructed node pair (i,j), i
%   should correspond to a node for rows and j corresponds to columns. 
%


enbox = htree.enbox;
row_children = row_htree.children;
col_children = col_htree.children;
root = length(htree.parent);

[NFList, FFList] = H2_interaction(root, root);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [near, far] = H2_interaction(p1, p2)
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
if (isnan(row_children(p1, 1)) && isnan(col_children(p2, 1)))
    near = [p1, p2];
    far = [];
    return ;
end

%   p1 leaf and p2 nonleaf
if (isnan(row_children(p1, 1)) && ~isnan(col_children(p2, 1)))
    nchild = find(~isnan(col_children(p2,:)), 1, 'last');
    for i = 1 : nchild
        [nn, ff] = H2_interaction(p1, col_children(p2, i));
        near = vertcat(near, nn);
        far  = vertcat(far , ff);
    end
    return ;
end

%   p1 nonleaf and p2 leaf
if (~isnan(row_children(p1, 1)) && isnan(col_children(p2,1)))
    nchild = find(~isnan(row_children(p1,:)), 1, 'last');
    for i = 1 : nchild
        [nn, ff] = H2_interaction(row_children(p1, i), p2);
        near = vertcat(near, nn);
        far  = vertcat(far, ff);
    end
    return ;
end

%   both p1 and p2 nonleaf 
if (~isnan(row_children(p1, 1)) && ~isnan(col_children(p2,1)))
    nchild1 = find(~isnan(row_children(p1,:)), 1, 'last');
    nchild2 = find(~isnan(col_children(p2,:)), 1, 'last');
    for i = 1 : nchild1
        for j = 1 : nchild2
            [nn, ff] = H2_interaction(row_children(p1, i), col_children(p2, j));
            near = vertcat(near, nn);
            far  = vertcat(far, ff);
        end
    end
end

end
end