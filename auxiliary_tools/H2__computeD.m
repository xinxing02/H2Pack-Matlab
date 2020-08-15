function D = H2__computeD(h2mat, htree, node1, node2)

    %   Basic info of tree
    cluster  = htree.cluster;
    coord    = htree.coord;
    kernel   = h2mat.kernel;
    idx1 = cluster(node1, 1): cluster(node1, 2);
    idx2 = cluster(node2, 1): cluster(node2, 2);    
    
    if node1 == node2
        D = kernel(coord(idx1,:));
    else
        D = kernel({coord(idx1, :), coord(idx2, :)});
    end
end