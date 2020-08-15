function A = H2__denseform(h2mat, htree)
%
%   Construct the dense matrix form of the H2 representation
%

    %%  Basic information of h2mat
    near = h2mat.near;
    far = h2mat.far;
    U = H2__explicitU(h2mat, htree);
    kdim = h2mat.kdim;

    %%  Basic information of htree
    leafnode = htree.leafnode;
    nodelvl  = htree.nodelvl;
    root     = htree.root;
    cluster  = htree.cluster;
    mcluster = [kdim * cluster(:,1) - kdim + 1, kdim * cluster(:, 2)];
    
    %%  Initialization
    N = mcluster(root,2) - mcluster(root,1);
    A = zeros(N);
    
    if h2mat.JIT == false
        %%  Diagonal blocks
        for i = 1  : length(leafnode)
            idx = mcluster(leafnode(i),1) : mcluster(leafnode(i),2);
            A(idx, idx) = h2mat.D{leafnode(i), leafnode(i)};
        end

        %%  Dense in-admissible blocks
        for i = 1 : size(near, 1)
            c1 = near(i,1);
            c2 = near(i,2);
            idx1 = mcluster(c1, 1): mcluster(c1, 2);
            idx2 = mcluster(c2, 1): mcluster(c2, 2); 
            A(idx1, idx2) = h2mat.D{c1,c2};
            A(idx2, idx1) = h2mat.D{c1,c2}';
        end

        %%  Dense admissible blocks
        for i = 1 : size(far, 1)
            c1 = far(i,1);
            c2 = far(i,2);
            idx1 = mcluster(c1, 1): mcluster(c1, 2);
            idx2 = mcluster(c2, 1): mcluster(c2, 2);    
            if nodelvl(c1) == nodelvl(c2)
                %   compression on both sides
                A(idx1, idx2) = U{c1} * h2mat.B{c1,c2}  * U{c2}';            
            elseif nodelvl(c1) > nodelvl(c2)  
                %   c2 is the leafnode at higher level. only compress on c1's side
                A(idx1, idx2) = U{c1} * h2mat.B{c1,c2};
            else
                %   c1 is the leafnode at higher level. only compress on c2's side
                A(idx1, idx2) = h2mat.B{c1,c2} * U{c2}';
            end
            A(idx2, idx1) = A(idx1, idx2)';
        end
    else
        coord = htree.coord;
        kernel = h2mat.kernel;
        
        %%  Diagonal blocks
        for i = 1  : length(leafnode)
            idx = cluster(leafnode(i), 1) : cluster(leafnode(i), 2);
            midx = mcluster(leafnode(i),1) : mcluster(leafnode(i),2);
            A(midx, midx) = kernel(coord(idx, :));
        end

        %%  Dense in-admissible blocks
        for i = 1 : size(near, 1)
            c1 = near(i,1);
            c2 = near(i,2);
            idx1 = cluster(c1, 1): cluster(c1, 2);
            idx2 = cluster(c2, 1): cluster(c2, 2);
            midx1 = mcluster(c1, 1): mcluster(c1, 2);
            midx2 = mcluster(c2, 1): mcluster(c2, 2);
            A(midx1, midx2) = kernel({coord(idx1, :), coord(idx2, :)});
            A(midx2, midx1) = A(midx1, midx2)';            
        end

        %%  Dense admissible blocks
        for i = 1 : size(far, 1)
            c1 = far(i,1);
            c2 = far(i,2);
            idx1 = cluster(c1, 1): cluster(c1, 2);
            idx2 = cluster(c2, 1): cluster(c2, 2);
            midx1 = mcluster(c1, 1): mcluster(c1, 2);
            midx2 = mcluster(c2, 1): mcluster(c2, 2);   
            if nodelvl(c1) == nodelvl(c2)
                %   compression on both sides
                A(midx1, midx2) = U{c1} * kernel({coord(h2mat.I{c1},:), coord(h2mat.I{c2},:)}) * U{c2}';            
            elseif nodelvl(c1) > nodelvl(c2)  
                %   c2 is the leafnode at higher level. only compress on c1's side
                A(midx1, midx2) = U{c1} * kernel({coord(h2mat.I{c1},:), coord(idx2, :)});
            else
                %   c1 is the leafnode at higher level. only compress on c2's side
                A(midx1, midx2) = kernel({coord(idx1, :), coord(h2mat.I{c2}, :)}) * U{c2}';
            end
            A(midx2, midx1) = A(midx1, midx2)';
        end
    end
end
        
    