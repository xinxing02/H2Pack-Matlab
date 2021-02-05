function u = H2_matvec(h2mat, htree, vec)
%
%   Compute A * u with H2 representation (D,U,B) for A. 
%
    
    %%  Basic information of the h2mat
    near = h2mat.near;
    far  = h2mat.far;
    kdim = h2mat.kdim;
    minlvl = h2mat.minlvl;
    kernel = h2mat.kernel;

    %%  Permute the input vector 
    %vec(htree.permutation, :) = vec;
    vec0 = vec;
    npt = length(htree.permutation);
    for i = 1 : npt
        j     = htree.permutation(i);
        sidx0 = (i - 1) * kdim + 1;
        eidx0 =  i      * kdim;
        sidx1 = (j - 1) * kdim + 1;
        eidx1 =  j      * kdim;
        vec(sidx1 : eidx1, :) = vec0(sidx0 : eidx0, :);
    end
    
    %%   Basic information of the htree
    children = htree.children;
    level    = htree.level;
    nlevel   = length(level);
    nodelvl  = htree.nodelvl;    
    leafnode = htree.leafnode;   
    coord    = htree.coord;
    cluster  = htree.cluster;
    mcluster = [kdim*cluster(:,1) - (kdim-1), kdim*cluster(:,2)];
    nnode    = htree.nnode;
    
    %%  Initialization
    nvec = size(vec, 2);
    u = zeros(size(vec));   
    rk = zeros(nnode, 1);
    
    %%   Upward Sweep (calculate the U_i' * y_i    
    y = cell(size(mcluster, 1), 1);
    for i = nlevel : -1 : minlvl
        for j = 1 : length(level{i})
            node = level{i}(j);
            child_node = children(node, ~isnan(children(node, :)));
            if isempty(child_node)
                y{node} = h2mat.U{node}' * vec(mcluster(node,1):mcluster(node,2), :);
            else
                y{node} = h2mat.U{node}' * vertcat(y{child_node});
            end
            rk(node) = size(h2mat.U{node},2);
        end
    end
    
    %%  Intermediate Sweep (calculate B_ij * U_j' * y_j)
    inter_u = cell(size(mcluster, 1),1);
    if h2mat.JIT == false        
        %  B precomputed
        for i = 1 : size(far, 1)
            c1 = far(i,1);
            c2 = far(i,2);

            if size(inter_u{c1}, 2) == 0
                inter_u{c1} = zeros(size(h2mat.U{c1}, 2), nvec);
            end
            if size(inter_u{c2}, 2) == 0
                inter_u{c2} = zeros(size(h2mat.U{c2}, 2), nvec);
            end

            if nodelvl(c1) == nodelvl(c2)
                %   compression on both sides
                inter_u{c1} = inter_u{c1} + h2mat.B{c1,c2} * y{c2};
                inter_u{c2} = inter_u{c2} + h2mat.B{c1,c2}'* y{c1};
            elseif nodelvl(c1) > nodelvl(c2)  
                %   c2 is the leafnode at higher level. only compress on c1's side
                inter_u{c1} = inter_u{c1} + h2mat.B{c1,c2} * vec(mcluster(c2,1):mcluster(c2,2), :);         
                u(mcluster(c2,1):mcluster(c2,2), :) = u(mcluster(c2,1):mcluster(c2,2), :) + ...
                    h2mat.B{c1,c2}' * y{c1}; 
            else
                %   c1 is the leafnode at higher level. only compress on c2's side                    
                u(mcluster(c1,1):mcluster(c1,2), :) = u(mcluster(c1,1):mcluster(c1,2), :) + ...
                    h2mat.B{c1,c2} * y{c2};
                inter_u{c2} = inter_u{c2} + h2mat.B{c1,c2}' * vec(mcluster(c1,1):mcluster(c1,2), :);
            end
        end
    else
        %  B dynamically computed        
        for i = 1 : size(far, 1)
            c1 = far(i,1);
            c2 = far(i,2);

            %   calculate Bij
            idx1 = cluster(c1, 1): cluster(c1, 2);
            idx2 = cluster(c2, 1): cluster(c2, 2);    
            if nodelvl(c1) == nodelvl(c2)
                %   compression on both sides
                tmpB = kernel({coord(h2mat.I{c1},:), coord(h2mat.I{c2},:)}); 
            elseif nodelvl(c1) > nodelvl(c2)  
                %   c2 is the leafnode at higher level. only compress on c1's side
                tmpB = kernel({coord(h2mat.I{c1},:), coord(idx2, :)});
            else
                %   c1 is the leafnode at higher level. only compress on c2's side
                tmpB = kernel({coord(idx1, :), coord(h2mat.I{c2}, :)});
            end

            %   initialize intermediate vector
            if size(inter_u{c1}, 2) == 0
                inter_u{c1} = zeros(size(h2mat.U{c1}, 2), nvec);
            end
            if size(inter_u{c2}, 2) == 0
                inter_u{c2} = zeros(size(h2mat.U{c2}, 2), nvec);
            end

            %   horizontal matvec        
            if nodelvl(c1) == nodelvl(c2)
                %   compression on both sides
                inter_u{c1} = inter_u{c1} + tmpB * y{c2};
                inter_u{c2} = inter_u{c2} + tmpB'* y{c1};
            elseif nodelvl(c1) > nodelvl(c2)  
                %   c2 is the leafnode at higher level. only compress on c1's side
                inter_u{c1} = inter_u{c1} + tmpB * vec(mcluster(c2,1):mcluster(c2,2), :);         
                u(mcluster(c2,1):mcluster(c2,2), :) = u(mcluster(c2,1):mcluster(c2,2), :) + ...
                    tmpB' * y{c1}; 
            else
                %   c1 is the leafnode at higher level. only compress on c2's side                    
                u(mcluster(c1,1):mcluster(c1,2), :) = u(mcluster(c1,1):mcluster(c1,2), :) + ...
                    tmpB * y{c2};
                inter_u{c2} = inter_u{c2} + tmpB' * vec(mcluster(c1,1):mcluster(c1,2), :);
            end
        end
    end
        
    
    %%  Down Sweep (calculate U_i * B_ij * U_j' * y_j)
    for i = minlvl : nlevel
        for j = 1 : length(level{i})
            node = level{i}(j); 
            if isempty(inter_u{node}); continue; end
            inter_u{node} = h2mat.U{node} * inter_u{node}; 
            child_node = children(node, ~isnan(children(node, :)));                                                      
            if isempty(child_node)  %   leaf node
                u(mcluster(node,1):mcluster(node,2), :) = u(mcluster(node,1):mcluster(node,2), :) + inter_u{node};
            else      %    non-leaf node      
                nrows = rk(child_node)';
%                 nrows = arrayfun(@(k)size(h2mat.U{k},2), child_node);
                offrow = [1 cumsum(nrows)+1];
                for k = 1 : length(child_node)
                    if isempty(inter_u{child_node(k)})
                        inter_u{child_node(k)} = inter_u{node}(offrow(k):offrow(k+1)-1, :);
                    else
                        inter_u{child_node(k)} = inter_u{child_node(k), :} + inter_u{node}(offrow(k):offrow(k+1)-1, :);
                    end
                end
            end
        end
    end 
    
    %%   Dense blocks
    if h2mat.JIT == false
        %   D precomputed
        for i = 1 : length(leafnode)
            node = leafnode(i);
            idx = mcluster(node, 1) : mcluster(node, 2);
            u(idx, :) = u(idx, :) + h2mat.D{node, node} * vec(idx, :);
        end

        for i = 1 : size(near, 1)
            c1 = near(i, 1);
            c2 = near(i, 2);
            idx1 = mcluster(c1, 1):mcluster(c1,2);
            idx2 = mcluster(c2, 1):mcluster(c2,2);
            u(idx1, :) = u(idx1, :) + h2mat.D{c1, c2} * vec(idx2, :);
            u(idx2, :) = u(idx2, :) + h2mat.D{c1, c2}' * vec(idx1, :);
        end
    else
        %   D dynamically computed
        for i = 1 : length(leafnode)
            node = leafnode(i);
            idx = cluster(node, 1) : cluster(node, 2);
            midx = mcluster(node, 1) : mcluster(node, 2);
            tmpD = kernel(coord(idx, :));
            u(midx, :) = u(midx, :) + tmpD * vec(midx, :);
        end
  
        for i = 1 : size(near, 1)
            c1 = near(i, 1);
            c2 = near(i, 2);
            idx1 = cluster(c1, 1): cluster(c1, 2);
            idx2 = cluster(c2, 1): cluster(c2, 2);
            tmpD = kernel({coord(idx1, :), coord(idx2, :)});

            midx1 = mcluster(c1, 1):mcluster(c1,2);
            midx2 = mcluster(c2, 1):mcluster(c2,2);
            u(midx1, :) = u(midx1, :) + tmpD * vec(midx2, :);
            u(midx2, :) = u(midx2, :) + tmpD' * vec(midx1, :);
        end
    end
    
    %%  Permute the output vector 
    %u = u(htree.permutation, :);
    u0 = u;
    npt = length(htree.permutation);
    for i = 1 : npt
        j     = htree.permutation(i);
        sidx0 = (i - 1) * kdim + 1;
        eidx0 =  i      * kdim;
        sidx1 = (j - 1) * kdim + 1;
        eidx1 =  j      * kdim;
        u(sidx0 : eidx0, :) = u0(sidx1 : eidx1, :);
    end
end
