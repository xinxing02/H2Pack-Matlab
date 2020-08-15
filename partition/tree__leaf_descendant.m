function nodelist = tree__leaf_descendant(htree)

    children = htree.children;
    level = htree.level;
    nlevel = length(level);
    parent = htree.parent;
    
    nodelist = cell(length(parent), 1);
    
    for i = nlevel : -1 : 1
        for j = 1 : length(level{i})
            node = level{i}(j);
            child_node = children(node, ~isnan(children(node, :)));
            
            if isempty(child_node)
                nodelist{node} = [node];
            end
            
            if ~isempty(child_node)
                nodelist{node} = horzcat(nodelist{child_node});
            end
        end
    end
        