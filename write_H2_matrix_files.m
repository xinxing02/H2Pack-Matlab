function write_H2_matrix_files(htree, h2mat, metadata_fname, binary_fname)
    metadata_fid = fopen(metadata_fname, 'w');
    binary_fid   = fopen(binary_fname,   'w');

    %% 1. Metadata: H2 / HSS common part
    n_point = size(htree.coord, 1);
    if (h2mat.alpha > 0.99)
        is_HSS = 0;
    else
        is_HSS = 1;
    end
    n_near = size(h2mat.near, 1);
    n_far  = size(h2mat.far,  1);
    pt_dim = size(htree.coord, 2);
    max_child = 2^pt_dim;
    kmat_size = h2mat.kdim * n_point;
    fprintf(metadata_fid, '%d\n', pt_dim);          % C.1 dim_point
    fprintf(metadata_fid, '%d\n', h2mat.kdim);      % C.2 dim_kernel
    fprintf(metadata_fid, '%d\n', n_point);         % C.3 num_point
    fprintf(metadata_fid, '%d\n', kmat_size);       % A.1 nrow_matrix
    fprintf(metadata_fid, '%d\n', kmat_size);       % A.2 ncol_matrix
    fprintf(metadata_fid, '%d\n', 1);               % A.3 is_symmetric
    fprintf(metadata_fid, '%d\n', htree.nnode);     % A.4 num_node_row
    fprintf(metadata_fid, '%d\n', htree.nnode);     % A.5 num_node_col
    fprintf(metadata_fid, '%d\n', htree.root-1);    % A.6 root_node_row
    fprintf(metadata_fid, '%d\n', htree.root-1);    % A.7 root_node_col
    fprintf(metadata_fid, '%d\n', htree.nlevel);    % A.8 num_level_row
    fprintf(metadata_fid, '%d\n', htree.nlevel);    % A.9 num_level_col
    fprintf(metadata_fid, '%d\n', is_HSS);          % C.4 is_HSS
    fprintf(metadata_fid, '%d\n', h2mat.minlvl-1);  % C.5 min_adm_level
    fprintf(metadata_fid, '%d\n', n_near);          % A.14 num_inadmissible_blocks - n_leaf_node
    fprintf(metadata_fid, '%d\n', n_far);           % A.15 num_admissible_blocks
    has_part_adm = 0;
    for i = 1 : n_far
        node0  = h2mat.far(i, 1);
        node1  = h2mat.far(i, 2);
        level0 = htree.nodelvl(node0);
        level1 = htree.nodelvl(node1);
        if (level0 ~= level1)
            has_part_adm = 1;
            break;
        end
    end
    fprintf(metadata_fid, '%d\n', has_part_adm);    % A.16 has_partial_adm_blocks
    cluster = htree.cluster;
    nodelvl = htree.nodelvl;

    %% 2. Metadata: partitioning tree
    % A.10 nodes_row; A.11 nodes_col == NULL since H2 matrix is symmetric
    for i = 1 : htree.nnode
        node_children = htree.children(i, :);
        n_child = sum(~isnan(node_children));
        fprintf(metadata_fid, '%6d ', i - 1);               % A.10.1 index
        fprintf(metadata_fid, '%2d ', nodelvl(i) - 1);      % A.10.2 level
        fprintf(metadata_fid, '%8d ', cluster(i, 1) - 1);   % A.10.3 cluster_head
        fprintf(metadata_fid, '%8d ', cluster(i, 2) - 1);   % A.10.4 cluster_tail
        fprintf(metadata_fid, '%2d ', n_child);             % A.10.5 num_children
        % A.10.6 children
        for j = 1 : n_child
            fprintf(metadata_fid, '%8d ', node_children(j) - 1);
        end
        for j = n_child + 1 : max_child
            fprintf(metadata_fid, '-1 ');
        end
        fprintf(metadata_fid, '\n');
    end

    %% 3. Metadata & binary data: U matrices
    %% A.12 basis_matrices_row (A.13 ignored since H2 matrix is symmetric)
    for i = 1 : htree.nnode
        tmpU = h2mat.U{i};
        [U_nrow, U_ncol] = size(tmpU);
        fprintf(metadata_fid, '%6d ',  i - 1);  % A.12.1 node
        fprintf(metadata_fid, '%5d ',  U_nrow); % A.12.2 num_row
        fprintf(metadata_fid, '%5d\n', U_ncol); % A.12.3 num_col
        fwrite(binary_fid, tmpU', 'double');    % B.1; B.2 == NULL since H2 matrix is symmetric
    end

    %% 4. Metadata & binary data: B matrices
    for i = 1 : n_far
        node0  = h2mat.far(i, 1);
        node1  = h2mat.far(i, 2);
        level0 = htree.nodelvl(node0);
        level1 = htree.nodelvl(node0);
        if (level0 == level1)
            is_part_adm = 0;
        else
            is_part_adm = 1;
        end
        if (h2mat.JIT == 1)
            idx0 = htree.cluster(node0, 1) : htree.cluster(node0, 2);
            idx1 = htree.cluster(node1, 1) : htree.cluster(node1, 2);  
            if htree.nodelvl(node0) == htree.nodelvl(node1)
                tmpB = h2mat.kernel({htree.coord(h2mat.I{node0},:), htree.coord(h2mat.I{node1},:)}); 
            elseif htree.nodelvl(node0) > htree.nodelvl(node1)  
                tmpB = h2mat.kernel({htree.coord(h2mat.I{node0},:), htree.coord(idx1, :)});
            else
                tmpB = h2mat.kernel({htree.coord(idx0, :), htree.coord(h2mat.I{node1}, :)});
            end
        else
            tmpB = h2mat.B{node0, node1};
        end
        [B_nrow, B_ncol] = size(tmpB);
        fprintf(metadata_fid, '%6d ', node0 - 1);   % A.17.1 node_row
        fprintf(metadata_fid, '%6d ', node1 - 1);   % A.17.2 node_col
        fprintf(metadata_fid, '%5d ', B_nrow);      % A.17.3 num_row
        fprintf(metadata_fid, '%5d ', B_ncol);      % A.17.4 num_col
        fprintf(metadata_fid, '%d\n', is_part_adm); % A.17.5 is_part_adm
        fwrite(binary_fid, tmpB', 'double');        % B.3
    end

    %% 5. Metadata & binary data: D matrices
    for i = 1 : length(htree.leafnode)
        node = htree.leafnode(i);
        if (h2mat.JIT == 1)
            idx  = htree.cluster(node, 1) : htree.cluster(node, 2);
            tmpD = h2mat.kernel(htree.coord(idx, :));
        else
            tmpD = h2mat.D{node, node};
        end
        [D_nrow, D_ncol] = size(tmpD);
        fprintf(metadata_fid, '%6d ',  node - 1);   % A.18.1 node_row
        fprintf(metadata_fid, '%6d ',  node - 1);   % A.18.2 node_col
        fprintf(metadata_fid, '%5d ',  D_nrow);     % A.18.3 num_row
        fprintf(metadata_fid, '%5d\n', D_ncol);     % A.18.4 num_col
        fwrite(binary_fid, tmpD', 'double');        % B.4
    end
    for i = 1 : n_near
        node0 = h2mat.near(i, 1);
        node1 = h2mat.near(i, 2);
        if (h2mat.JIT == 1)
            idx0 = htree.cluster(node0, 1) : htree.cluster(node0, 2);
            idx1 = htree.cluster(node1, 1) : htree.cluster(node1, 2);
            tmpD = h2mat.kernel({htree.coord(idx0, :), htree.coord(idx1, :)});
        else
            tmpD = h2mat.D{node0, node1};
        end
        [D_nrow, D_ncol] = size(tmpD);
        fprintf(metadata_fid, '%6d ',  node0 - 1);  % A.18.1 node_row
        fprintf(metadata_fid, '%6d ',  node1 - 1);  % A.18.2 node_col
        fprintf(metadata_fid, '%5d ',  D_nrow);     % A.18.3 num_row
        fprintf(metadata_fid, '%5d\n', D_ncol);     % A.18.4 num_col
        fwrite(binary_fid, tmpD', 'double');        % B.4
    end

    %% 6. Other necessary information for H2Pack
    fprintf(metadata_fid, '%d\n', htree.minsize);   % C.6 max_leaf_points
    fprintf(metadata_fid, '%e\n', h2mat.par(1));    % C.7 QR_stop_tol
    fprintf(metadata_fid, '%d\n', 1);               % C.8 has_skeleton_points
    % The permutation array in H2Pack-Matlab is the inverse function of H2Pack permutation array
    coord0 = zeros(n_point, pt_dim);
    perm   = htree.permutation;
    perm1  = zeros(size(perm));
    for i = 1 : n_point
        coord0(i, :) = htree.coord(perm(i), :);
        perm1(perm(i)) = i - 1;
    end
    % C.9 point_coordinate
    % Cast it to uint64_t so reading it from text file in MATLAB won't be so slow
    for i = 1 : n_point
        coord_i = typecast(coord0(i, :), 'uint64');
        for j = 1 : pt_dim
            fprintf(metadata_fid, '%x ', coord_i(j));
        end
        fprintf(metadata_fid, '\n');
    end
    % C.10 permutation_array
    for i = 1 : n_point
        fprintf(metadata_fid, '%d\n', perm1(i));
    end
    % C.11 skeleton_point
    for i = 1 : htree.nnode
        tmpI = h2mat.I{i};
        fprintf(metadata_fid, '%6d ', i - 1);           % C.11.1 node
        fprintf(metadata_fid, '%6d ', length(tmpI));    % C.11.2 num_skeleton_point_indices
        % C.11.3 skeleton_point_indices
        for j = 1 : length(tmpI)
            fprintf(metadata_fid, '%d ', tmpI(j) - 1);
        end
        fprintf(metadata_fid, '\n');
    end

    fclose(metadata_fid);
    fclose(binary_fid);
end