function u = H2_matvec_sub(h2mat, htree, vec, rowidx, colidx)
%
%   Compute y = A(rowidx, colidx) * vec with an H2 matrix A
%    
    %%  Basic information of the h2mat
    near = h2mat.near;
    far  = h2mat.far;
    kdim = h2mat.kdim;
    minlvl = h2mat.minlvl;
    kernel = h2mat.kernel;
    
    %%   Basic information of the htree
    parent   = htree.parent;
    children = htree.children;
    level    = htree.level;
    nlevel   = length(level);
    nodelvl  = htree.nodelvl;    
    leafnode = htree.leafnode;   
    coord    = htree.coord;
    cluster  = htree.cluster;
    mcluster = [kdim*cluster(:,1) - (kdim-1), kdim*cluster(:,2)];
    nnode    = htree.nnode;
    htree.mcluster = mcluster;
     
    %%  Input vector check
    nvec = size(vec, 2);
    if size(vec, 1) ~= length(colidx)
        fprintf("Input vector length is inconsistent with the provided row indices\n");
    end
    [colidx, ii] = sort(colidx, 'ascend');
    vec = vec(ii, :);        
    
    %%  Locate leafnodes that contain rowidx and colidx
    row_leafnode = tree__idx2leafnode(htree, rowidx);
    col_leafnode = tree__idx2leafnode(htree, colidx);
    row_htree = tree__prune(htree, row_leafnode);
    col_htree = tree__prune(htree, col_leafnode);    
    
    %%  Construct interaction list for the interaction between the two sets of leaf nodes
    [NFList, FFList] = H2__interaction_list(htree, row_htree, col_htree, h2mat.alpha);
    
    %%  Properly expand the vector indexed by row_idx to the larger vector indexed by row_leafnode
    vec_in = cell(nnode, 1);
    for i = 1 : length(col_leafnode)
        node = col_leafnode(i);
        idx = find( colidx < mcluster(node, 2) + 1 & colidx > mcluster(node, 1) -1 );
        vec_in{node} = zeros(mcluster(node, 2) - mcluster(node, 1) + 1, nvec);
        vec_in{node}( colidx(idx) - mcluster(node, 1) +1 , :) = vec(idx, :);
    end
            
    %%  Initialize output vectors
    vec_out = cell(nnode, 1);
    for i = 1 : length(row_leafnode)
        node = row_leafnode(i);
        vec_out{node} = zeros(mcluster(node, 2) - mcluster(node, 1) + 1, nvec);
    end
    

    
    %%   Upward Sweep (calculate the U_i' * y_i    
    upward_minlvl = min(nodelvl(FFList(:, 2)));
    col_level     = col_htree.level;
    col_children  = col_htree.children;
    col_nlevel    = col_htree.nlevel;
    
    y = cell(nnode, 1);
    y_size = zeros(1, nnode);   %   number of rows for y at node i
    for i = 1 : nnode
        y_size(i) = size(h2mat.U{i},2);
    end
    
    for i = col_nlevel : -1 : upward_minlvl
        for j = 1 : length(col_level{i})
            node = col_level{i}(j);
            col_childnode = col_children(node, ~isnan(col_children(node, :)));
            if isempty(col_childnode)                
                y{node} = h2mat.U{node}' * vec_in{node};
            else
                all_childnode = children(node, ~isnan(children(node, :)));
                %   note: this assumed ascending order of clusters defined
                %   by children nodes in children(i, :). 
                offset = [1, cumsum(y_size(all_childnode))+1];
                idx = cell(length(col_childnode), 1);
                for k = 1 : length(col_childnode)                    
                    cnode = col_childnode(k);
                    ii = find(all_childnode == cnode, 1, 'first');
                    idx{k} = offset(ii) : (offset(ii+1)-1);                    
                end
                idx = horzcat(idx{:});
                y{node} = h2mat.U{node}(idx, :)' * vertcat(y{col_childnode});
            end
        end
    end
    
    %%  Intermediate Sweep (calculate B_ij * U_j' * y_j)    
    inter_u = cell(nnode, nvec);
    for i = 1 : size(FFList, 1)
        c1 = FFList(i,1);
        c2 = FFList(i,2);
        
        if h2mat.JIT == false
            if c1 > c2
                tmpB = h2mat.B{c2,c1}';
            else
                tmpB = h2mat.B{c1,c2};
            end
        else
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
        end
        
        if size(inter_u{c1}, 2) == 0
            inter_u{c1} = zeros(size(h2mat.U{c1}, 2), nvec);
        end
        
        if nodelvl(c1) == nodelvl(c2)
            %   compression on both sides
            inter_u{c1} = inter_u{c1} + tmpB * y{c2};
        elseif nodelvl(c1) > nodelvl(c2)  
            %   c2 is the leafnode at higher level. only compress on c1's side
            inter_u{c1} = inter_u{c1} + tmpB * vec_in{c2};         
        else
            %   c1 is the leafnode at higher level. only compress on c2's side                    
            vec_out{c1} = vec_out{c1} + tmpB * y{c2};
        end
    end      
    
    %%  Down Sweep (calculate U_i * B_ij * U_j' * y_j)
    downward_minlvl = min(nodelvl(FFList(:, 1)));
    row_level     = row_htree.level;
    row_children  = row_htree.children;
    row_nlevel    = row_htree.nlevel;
    for i = downward_minlvl : row_nlevel
        for j = 1 : length(row_level{i})
            node = row_level{i}(j); 
            if size(inter_u{node},2)==0; continue; end            
            
            row_childnode = row_children(node, ~isnan(row_children(node, :)));                                                                  
            if isempty(row_childnode)
                vec_out{node} = vec_out{node} + h2mat.U{node} * inter_u{node}; 
            else
                tmpU = h2mat.U{node};
                all_childnode = children(node, ~isnan(children(node, :)));
                %   note: this assumed ascending order of clusters defined
                %   by children nodes in children(i, :). 
                offset = [1, cumsum(y_size(all_childnode))+1];
                for k = 1 : length(row_childnode)
                    cnode = row_childnode(k);
                    ii = find(all_childnode == cnode, 1, 'first');
                    idx = offset(ii) : (offset(ii+1)-1);
                    if size(inter_u{cnode}, 1) == 0
                        inter_u{cnode} = tmpU(idx, :) * inter_u{node};
                    else
                        inter_u{cnode} = inter_u{cnode} + tmpU(idx, :) * inter_u{node};
                    end
                end
            end
        end
    end 
    
    %%  Dense block matvec               
    for i = 1 : size(NFList, 1)
        c1 = NFList(i, 1);
        c2 = NFList(i, 2);
        
        if h2mat.JIT == false
            if c1 > c2
                tmpD = h2mat.D{c2, c1}';
            else
                tmpD = h2mat.D{c1, c2};
            end
        else
            idx1 = cluster(c1, 1): cluster(c1, 2);
            idx2 = cluster(c2, 1): cluster(c2, 2);
            if c1 == c2
                tmpD = kernel(coord(idx1, :));
            else
                tmpD = kernel({coord(idx1, :), coord(idx2, :)});
            end
        end
        vec_out{c1} = vec_out{c1} + tmpD * vec_in{c2};
    end
   
    %%  Restrict vec_out from row_leafnode to row_idx
    [rowidx, perm] = sort(rowidx, 'ascend');
    u = zeros(length(rowidx), nvec);    
    for i = 1 : length(row_leafnode)
        node = row_leafnode(i);
        idx = find( rowidx < mcluster(node, 2) + 1 & rowidx > mcluster(node, 1) -1 );
        u(idx, :) = vec_out{node}( rowidx(idx) - mcluster(node, 1) + 1, :);
    end
    u(perm, :) = u;
end
