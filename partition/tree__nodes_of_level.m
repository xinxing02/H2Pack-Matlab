function lvl = tree__nodes_of_level(parent)
%
%   Parse the parent vector to get nodes' indices at each level of the tree.
%   
    ntree = length(parent);
    flag = ones(ntree,1);
    lvl = {[]};
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
        if s > length(lvl)
            lvl{s} = [];
        end
        for j = 1 : s
            if flag(path(j)) == 0
                continue;
            end
            lvl{s+1-j} = [lvl{s+1-j}, path(j)];
            flag(path(j)) = 0;
        end
    end
    for i = 1 : length(lvl)
        lvl{i} = sort(lvl{i}, 'ascend');
    end
end
            
        