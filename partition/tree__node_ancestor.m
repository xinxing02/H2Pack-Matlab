function p = tree__node_ancestor(parent, node)
%
%   The list of the ancestor indices of 'node'; 
%
    if node == length(parent)
        p = [];       
    else
        p = [parent(node), tree__node_ancestor(parent, parent(node))];
    end      
end