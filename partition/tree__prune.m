function subtree = tree__prune(htree, nodes)
%
%   Given a partition tree 'htree' and a set of tree nodes, 'nodes', 
%   construct a subtree of htree that only contains 'nodes' and their
%   ancestors. 
%


parent   = htree.parent;
children = htree.children;
num_node = length(parent);
max_child = size(children, 2);

subtree_nodes = false(num_node, 1);
for i = 1 : length(nodes)
    node = nodes(i);
    ance = [node, tree__node_ancestor(parent, node)];
    subtree_nodes(ance) = true;    
end

parent(~subtree_nodes) = nan;
children(~subtree_nodes, :) = nan;
for i = 1 : num_node
    if subtree_nodes(i) == false
        continue;
    end
    childnode = children(i, ~isnan(children(i, :)));
    if isempty(childnode)
        continue;
    end
    childnode( ~subtree_nodes(childnode) ) = [];
    children(i, :) = [childnode, nan(1, max_child - length(childnode))];
end

level = htree.level;
nlevel = 0;
for i = 1 : length(level)
    level{i}(~subtree_nodes(level{i})) = [];
    if ~isempty(level{i})
        nlevel = i;
    end
end
nodelvl = htree.nodelvl;
nodelvl(~subtree_nodes) = nan;

subtree.parent   = parent;
subtree.children = children;
subtree.level    = level; 
subtree.nodelvl  = nodelvl;
subtree.leafnode = nodes;
subtree.root     = find(~(htree.parent>0));
subtree.nlevel   = nlevel;