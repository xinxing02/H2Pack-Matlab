function [htree, h2mat] = read_H2_matrix_files(metadata_fname, binary_fname)
    metadata_fid = fopen(metadata_fname);
    binary_fid   = fopen(binary_fname);

    %% 1. Metadata: H2 / HSS common part
    htree.dim       = fscanf(metadata_fid, '%d', 1);    % C.1 dim_point
    h2mat.kdim      = fscanf(metadata_fid, '%d', 1);    % C.2 dim_kernel
    htree.n_point   = fscanf(metadata_fid, '%d', 1);    % C.3 num_point
    htree.kmat_size = fscanf(metadata_fid, '%d', 1);    % A.1 nrow_matrix
    htree.kmat_size = fscanf(metadata_fid, '%d', 1);    % A.2 ncol_matrix
    htree.is_symm   = fscanf(metadata_fid, '%d', 1);    % A.3 is_symmetric
    htree.nnode     = fscanf(metadata_fid, '%d', 1);    % A.4 num_node_row
    htree.nnode     = fscanf(metadata_fid, '%d', 1);    % A.5 num_node_col
    htree.root      = fscanf(metadata_fid, '%d', 1)+1;  % A.6 root_node_row
    htree.root      = fscanf(metadata_fid, '%d', 1)+1;  % A.7 root_node_col
    htree.nlevel    = fscanf(metadata_fid, '%d', 1);    % A.8 num_level_row
    htree.nlevel    = fscanf(metadata_fid, '%d', 1);    % A.9 num_level_col
    h2mat.is_HSS    = fscanf(metadata_fid, '%d', 1);    % C.4 is_HSS
    h2mat.minlvl    = fscanf(metadata_fid, '%d', 1)+1;  % C.5 min_adm_level
    h2mat.n_near    = fscanf(metadata_fid, '%d', 1);    % A.14 num_inadmissible_blocks - n_leaf_node
    h2mat.n_far     = fscanf(metadata_fid, '%d', 1);    % A.15 num_admissible_blocks
    has_part_adm    = fscanf(metadata_fid, '%d', 1);    % A.16 has_partial_adm_blocks
    if (h2mat.is_HSS == 1)
        h2mat.alpha = 0;
    else
        h2mat.alpha = 1;
    end
    h2mat.JIT  = 0;
    h2mat.type = 'reltol';
    max_child  = 2^htree.dim;
    htree.parent   = zeros(htree.nnode, 1);
    htree.children = zeros(htree.nnode, max_child);
    htree.cluster  = zeros(htree.nnode, 2);
    htree.nodelvl  = zeros(htree.nnode, 1);
    htree.leafnode = zeros(htree.nnode, 1);
    htree.level    = cell(1, htree.nlevel);
    for i = 1 : htree.nlevel
        htree.level{i} = [];
    end
    n_leaf = 0;

    %% 2. Metadata: partitioning tree
    % A.10 nodes_row; A.11 nodes_col == NULL since H2 matrix is symmetric
    for i = 1 : htree.nnode
        node = fscanf(metadata_fid, '%d', 1) + 1;   % A.10.1 index
        lvl  = fscanf(metadata_fid, '%d', 1) + 1;   % A.10.2 level
        htree.nodelvl(node)    = lvl;
        htree.cluster(node, 1) = fscanf(metadata_fid, '%d', 1) + 1; % A.10.3 cluster_head
        htree.cluster(node, 2) = fscanf(metadata_fid, '%d', 1) + 1; % A.10.4 cluster_tail
        htree.level{lvl} = [htree.level{lvl} node];
        n_child = fscanf(metadata_fid, '%d', 1);    % A.10.5 num_children
        % A.10.6 children
        if (n_child == 0)
            n_leaf = n_leaf + 1;
            htree.leafnode(n_leaf) = node;
        end
        for j = 1 : max_child
            htree.children(node, j) = fscanf(metadata_fid, '%d', 1) + 1;
            if (htree.children(node, j) == 0)
                htree.children(node, j) = NaN; 
            else
                htree.parent(htree.children(node, j)) = node;
            end
        end
    end
    for i = 1 : htree.nlevel
        htree.level{i} = sort(htree.level{i}, 'ascend');
    end
    htree.n_leaf = n_leaf;
    htree.leafnode = htree.leafnode(1 : n_leaf);

    %% 3. Metadata & binary data: U matrices
    %% A.12 basis_matrices_row (A.13 ignored since H2 matrix is symmetric)
    h2mat.U = cell(htree.nnode, 1);
    for i = 1 : htree.nnode
        U_idx  = fscanf(metadata_fid, '%d', 1) + 1; % A.12.1 node
        U_nrow = fscanf(metadata_fid, '%d', 1);     % A.12.2 num_row
        U_ncol = fscanf(metadata_fid, '%d', 1);     % A.12.3 num_col
        tmpU   = fread(binary_fid, [U_ncol, U_nrow], 'double'); % B.1; B.2 == NULL since H2 matrix is symmetric
        h2mat.U{U_idx} = tmpU';
    end

    %% 4. Metadata & binary data: B matrices
    h2mat.far = zeros(h2mat.n_far, 2);
    h2mat.B   = cell(htree.nnode, htree.nnode);
    for i = 1 : h2mat.n_far
        node0       = fscanf(metadata_fid, '%d', 1) + 1;    % A.17.1 node_row
        node1       = fscanf(metadata_fid, '%d', 1) + 1;    % A.17.2 node_col
        B_nrow      = fscanf(metadata_fid, '%d', 1);        % A.17.3 num_row
        B_ncol      = fscanf(metadata_fid, '%d', 1);        % A.17.4 num_col
        is_part_adm = fscanf(metadata_fid, '%d', 1);        % A.17.5 is_part_adm
        tmpB        = fread(binary_fid, [B_ncol, B_nrow], 'double');    % B.3
        h2mat.far(i, 1) = node0;
        h2mat.far(i, 2) = node1;
        h2mat.B{node0, node1} = tmpB';
    end

    %% 5. Metadata & binary data: D matrices
    h2mat.near = zeros(h2mat.n_near, 2);
    h2mat.D    = cell(htree.nnode, htree.nnode);
    for i = 1 : htree.n_leaf
        node   = fscanf(metadata_fid, '%d', 1) + 1; % A.18.1 node_row
        node   = fscanf(metadata_fid, '%d', 1) + 1; % A.18.2 node_col
        D_nrow = fscanf(metadata_fid, '%d', 1);     % A.18.3 num_row
        D_ncol = fscanf(metadata_fid, '%d', 1);     % A.18.4 num_col
        tmpD   = fread(binary_fid, [D_ncol, D_nrow], 'double'); % B.4
        h2mat.D{node, node} = tmpD';
    end
    for i = 1 : h2mat.n_near
        node0  = fscanf(metadata_fid, '%d', 1) + 1; % A.18.1 node_row
        node1  = fscanf(metadata_fid, '%d', 1) + 1; % A.18.2 node_col
        D_nrow = fscanf(metadata_fid, '%d', 1);     % A.18.3 num_row
        D_ncol = fscanf(metadata_fid, '%d', 1);     % A.18.4 num_col
        tmpD   = fread(binary_fid, [D_ncol, D_nrow], 'double'); % B.4
        h2mat.near(i, 1) = node0;
        h2mat.near(i, 2) = node1;
        h2mat.D{node0, node1} = tmpD';
    end

    %% 6. Other necessary information for H2Pack
    htree.minsize = fscanf(metadata_fid, '%d', 1);  % C.6 max_leaf_points
    h2mat.par     = fscanf(metadata_fid, '%e', 1);  % C.7 QR_stop_tol
    has_skel      = fscanf(metadata_fid, '%d', 1);  % C.8 has_skeleton_points
    % C.9 point_coordinate
    % Cast it from uint64_t back to double
    coord0 = zeros(htree.n_point, htree.dim);
    for i = 1 : htree.n_point
        coord_i = zeros(1, htree.dim, 'uint64');
        for j = 1 : htree.dim
            coord_i(j) = fscanf(metadata_fid, '%lx', 1);
        end
        coord0(i, :) = typecast(coord_i, 'double');
    end
    prem0 = zeros(htree.n_point, 1);
    for i = 1 : htree.n_point
        prem0(i) = fscanf(metadata_fid, '%d', 1) + 1;
    end
    % The permutation array in H2Pack-Matlab is the inverse function of H2Pack permutation array
    perm  = zeros(htree.n_point, 1);
    coord = coord0;
    for i = 1 : htree.n_point
        perm(prem0(i)) = i;
        coord(i, :) = coord0(prem0(i), :);
    end
    htree.permutation = perm;
    htree.coord = coord;
    h2mat.I = cell(htree.nnode, 1);
    % C.11 skeleton_point
    if (has_skel == 1)
        for i = 1 : htree.nnode
            node   = fscanf(metadata_fid, '%d', 1) + 1;
            n_skel = fscanf(metadata_fid, '%d', 1);
            tmpI   = zeros(1, n_skel);
            for j = 1 : n_skel
                tmpI(j) = fscanf(metadata_fid, '%d', 1) + 1;
            end
            h2mat.I{node} = tmpI;
        end
    end

    fclose(metadata_fid);
    fclose(binary_fid);
end
