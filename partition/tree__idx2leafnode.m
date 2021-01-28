function node = tree__idx2leafnode(htree, idx)
%
%   Return the set of leafnodes that contain the given row/column indices
%
    leafnode = htree.leafnode;
    mcluster = htree.mcluster;
    mcluster_leaf = mcluster(leafnode, :);
    
    idx = sort(idx, 'ascend');
    flag = false(length(leafnode), 1);    
    for i = 1 : length(leafnode)
        flag(i) = ~isempty(find( (idx < mcluster_leaf(i,2)+1) & (idx > mcluster_leaf(i,1)-1), 1, 'first')); 
    end
    node = leafnode(flag); 
end