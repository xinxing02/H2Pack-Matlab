function storage = H2__storage_cost(h2mat, htree)   
    near = h2mat.near;
    far  = h2mat.far;
    kdim = h2mat.kdim;
    U    = h2mat.U;
    
    cluster = htree.cluster;
    mcluster = [kdim*cluster(:,1)-kdim+1, cluster(:,2)];
    leafnode = htree.leafnode;
    nodelvl  = htree.nodelvl;
    clstsize = mcluster(:,2) - mcluster(:,1) + 1;
    
    nD = 0; 
    nU = 0; 
    nB = 0;
    nadm = 0;

    %   inadmissible diagonal blocks
    for i = 1 : length(leafnode)        
        nD = nD + clstsize(leafnode(i))^2;
    end
    
    %   inadmissible off-diagonal blocks
    for i = 1 : size(near, 1)
        nD = nD + clstsize(near(i,1)) * clstsize(near(i,2));    
    end
    
    %   basis matrix U 
    for i = 1 : size(mcluster,1)
        nU = nU + numel(U{i});  
    end
        
    %   intermediate matrix B
    for i = 1 : size(far, 1)
        c1 = far(i,1);
        c2 = far(i,2);
        if nodelvl(c1) == nodelvl(c2)
            nB = nB + size(U{c1},2) * size(U{c2},2);
        elseif nodelvl(c1) < nodelvl(c2)
            % c2 is the leafnode at higher level. only compress on c1's side
            nB = nB + size(U{c1},2) * clstsize(c2);
        else
            % c1 is the leafnode at higher level. only compress on c2's side
            nB = nB + size(U{c2},2) * clstsize(c1);
        end
    end
    
    %   admissible blocks
    for i = 1 : size(far, 1)
        nadm = nadm + clstsize(far(i,1)) * clstsize(far(i,2));
    end  
    
    storage = [nD, nU, nB, nadm];
end    
    
        