function dimrank = H2__rank(h2mat, htree)
%
%   Return the rank of the associated H2 block row at each node
%   with form  [node, dim, rank];
%
    kdim = h2mat.kdim;
    U = h2mat.U;
    
    level = htree.level;
    nlevel = length(level);
    cluster = htree.cluster;
    mcluster = [kdim*cluster(:,1)-kdim+1, kdim*cluster(:,2)];
    
    dimrank = cell(nlevel, 1);
    for i = 1 : nlevel
        dimrank{i} = zeros(length(level{i}), 3);
        for j = 1 : length(level{i})
            node = level{i}(j);
            if isempty(U{node})
                dimrank{i}(j, :) = [node, mcluster(node,2)-mcluster(node,1)+1, inf];
            else
                dimrank{i}(j, :) = [node, mcluster(node,2)-mcluster(node,1)+1, size(U{node},2)];
            end
        end
    end
end
