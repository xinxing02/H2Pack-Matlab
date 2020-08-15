function h2mat = H2__denseblocks(h2mat, htree)

        
    %   Basic info of tree
    nodelvl  = htree.nodelvl;
    leafnode = htree.leafnode;
    cluster  = htree.cluster;
    coord    = htree.coord;
    kernel = h2mat.kernel;
    h2mat.D = [];
    h2mat.B = [];
    near = h2mat.near;
    far = h2mat.far;
    
    %   Initialization
    tnode = length(htree.nnode);
    D = cell(tnode);
    B = cell(tnode);



    %   diagonal blocks   
    for i = 1 : length(leafnode)
        node = leafnode(i);
        idx  = cluster(node,1) : cluster(node,2);
        D{node, node} = kernel(coord(idx,:));
    end

    %   off-diagonal in-admissible blocks
    for i = 1 : size(near, 1)
        c1 = near(i,1);
        c2 = near(i,2);
        idx1 = cluster(c1, 1): cluster(c1, 2);
        idx2 = cluster(c2, 1): cluster(c2, 2);
        D{c1, c2} = kernel({coord(idx1, :), coord(idx2, :)});
    end   

    %   intermediate blocks
    for i = 1 : size(far, 1)
        c1 = far(i,1);
        c2 = far(i,2);
        idx1 = cluster(c1, 1): cluster(c1, 2);
        idx2 = cluster(c2, 1): cluster(c2, 2);    
        if nodelvl(c1) == nodelvl(c2)
            %   compression on both sides
            B{c1,c2} = kernel({coord(h2mat.I{c1},:), coord(h2mat.I{c2},:)}); 
        elseif nodelvl(c1) > nodelvl(c2)  
            %   c2 is the leafnode at higher level. only compress on c1's side
            B{c1,c2} = kernel({coord(h2mat.I{c1},:), coord(idx2, :)});
        else
            %   c1 is the leafnode at higher level. only compress on c2's side
            B{c1,c2} = kernel({coord(idx1, :), coord(h2mat.I{c2}, :)});
        end
    end
    
    %
    h2mat.JIT = false;
    h2mat.D = D;
    h2mat.B = B;
end