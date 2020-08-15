function [diagblk, inadmblk, admblk] = H2__norm_block(A, htree, near, far)
    if isfield(htree, 'mcluster')
        mcluster = htree.mcluster;
    else
        mcluster = htree.cluster;
    end
    leafnode = htree.leafnode;
    diagblk = 0;
    inadmblk = 0;
    admblk = 0;
    
    %   diagonal block norm
    for i = 1 : size(leafnode, 1)
        node = leafnode(i);
        idx = mcluster(node,1):mcluster(node,2);
        diagblk = diagblk  + norm(A(idx, idx), 'fro')^2;
    end
    diagblk = sqrt(diagblk);
            
    %   inadmissible off-diagonal blocks
    for i = 1 : size(near, 1)
        nd1 = near(i,1);
        nd2 = near(i,2);
        idx1 = mcluster(nd1, 1):mcluster(nd1, 2);
        idx2 = mcluster(nd2, 1):mcluster(nd2, 2);
        inadmblk = inadmblk + 2 * norm(A(idx1, idx2), 'fro')^2;        
    end
    inadmblk = sqrt(inadmblk);
    
    %   admissible off-diagonal blocks
    for i = 1 : size(far, 1)
        nd1 = far(i,1);
        nd2 = far(i,2);
        idx1 = mcluster(nd1, 1):mcluster(nd1, 2);
        idx2 = mcluster(nd2, 1):mcluster(nd2, 2);
        admblk = admblk + 2 * norm(A(idx1, idx2), 'fro')^2;        
    end
    admblk = sqrt(admblk);
end
    
        