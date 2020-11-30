function htree = hierarchical_partition(coord, minSize, dim, box)

    %%   Check the enclosing box
    if nargin < 4
        box = zeros(2, dim);
        center = mean(coord, 1);
        rela_coord = abs(bsxfun(@minus, coord, center));
        L = max(rela_coord(:)) + 1e-10;
        box(1, :) = center - L * ones(1, dim);
        box(2, :) = 2 * L * ones(1, dim);
    end
    
    if min(abs(box(2,:))) < eps('double')
        disp 'point lying in a lower-dimension'
        htree = [];
        return ;
    end
    
    if size(box, 2) < dim
        disp 'Inconsistent box dimension with the given dim'
        htree = [];
        return ;
    end     

    %   minor modification to avoid rounding error
    box(1,:) = box(1,:) - 200*eps('double');
    box(2,:) = box(2,:) + 400*eps('double');

    [parent, children, cluster, enbox, permutation, ntree] = partition_bisec(coord, minSize, dim, box);
    htree.permutation = permutation;
    htree.parent   = parent;
    htree.children = children; 
    htree.cluster  = cluster;
    htree.minsize  = minSize;
    htree.level    = tree__nodes_of_level(parent);
    htree.nlevel   = length(htree.level);
    htree.nodelvl  = tree__level_of_node(parent);
    htree.leafnode = find(isnan(children(:,1)));
    htree.root     = find(~(parent>0));
    htree.enbox    = enbox;
    htree.nnode    = ntree;
    coord(permutation, :) = coord;
    htree.coord = coord;
end