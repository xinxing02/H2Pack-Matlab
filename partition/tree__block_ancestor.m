function pblock = tree__block_ancestor(htree, block)
    pblock = zeros(size(block));    
    for i = 1 : size(block, 1)
        c1 = block(i, 1);
        c2 = block(i, 2);
        %   c1 and c2 have to be different
        p1 = [flip(tree__node_ancestor(htree.parent, c1)), c1];
        p2 = [flip(tree__node_ancestor(htree.parent, c2)), c2];
        len = min(length(p1), length(p2));
        idx = find(p1(1:len) == p2(1:len), 1, 'last');
        pblock(i, :) = [p1(idx+1), p2(idx+1)];                
    end
end