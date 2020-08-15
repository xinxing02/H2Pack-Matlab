function exU = H2__explicitU(h2mat, htree)
%
%   Calculate the explicit U at each level of HSS/H2 form
%   
U = h2mat.U;
level    = htree.level;
children = htree.children;
nlevel   = length(level);
exU      = {[]};

for i = nlevel : -1 : 2
    for j = 1 : length(level{i})
        node = level{i}(j);
        child_node = children(node, ~isnan(children(node,:)));
        
        if (length(U) < node); continue; end
        if (isempty(U{node})); continue; end
        if isempty(child_node) %Leaf node
            exU{node} = U{node};
        else
            exU{node} = blkdiag(exU{child_node}) * U{node};
        end
    end
end