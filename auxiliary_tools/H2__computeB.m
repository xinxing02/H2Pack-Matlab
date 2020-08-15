function B = H2__computeB(h2mat, htree, node1, node2)
    %   Basic info of tree
    nodelvl  = htree.nodelvl;
    cluster  = htree.cluster;
    coord    = htree.coord;
    kernel   = h2mat.kernel;
    idx1 = cluster(node1, 1): cluster(node1, 2);
    idx2 = cluster(node2, 1): cluster(node2, 2);    
    if nodelvl(node1) == nodelvl(node2)
        %   compression on both sides
        B = kernel({coord(h2mat.I{node1},:), coord(h2mat.I{node2},:)}); 
    elseif nodelvl(node1) > nodelvl(node2)  
        %   node2 is the leafnode at higher level. only compress on node1's side
        B = kernel({coord(h2mat.I{node1},:), coord(idx2, :)});
    else
        %   node1 is the leafnode at higher level. only compress on node2's side
        B = kernel({coord(idx1, :), coord(h2mat.I{node2}, :)});
    end