function node = tree__level_of_node(parent)
%
%   Parse the parent vector to get the level index of each node.
%   
    ntree = length(parent);
    flag = ones(ntree,1);
    level = {[]};
    for i = 1 : ntree
        if flag(i) == 0
            continue;
        end
        path = [i];
        par = parent(i);
        while par > 0
            path = [path, par];
            par = parent(par);
        end
        s = length(path);
        if s > length(level)
            level{s} = [];
        end
        for j = 1 : s
            if flag(path(j)) == 0
                continue;
            end
            level{s+1-j} = [level{s+1-j}, path(j)];
            flag(path(j)) = 0;
        end
    end
    node = zeros(length(parent), 1);
    for i = 1 : length(level)
        node(level{i}) = i;
    end
end
            
        