function [htree, h2mat] = read_H2_matrix_files(metadata_fname, binary_fname)
    metadata_fid = fopen(metadata_fname);
    binary_fid   = fopen(binary_fname);

    %% 1. Metadata: H2 / HSS common part
    htree.dim       = fscanf(metadata_fid, '%d', 1);
    h2mat.kdim      = fscanf(metadata_fid, '%d', 1);
    htree.n_point   = fscanf(metadata_fid, '%d', 1);
    htree.kmat_size = fscanf(metadata_fid, '%d', 1);
    htree.is_symm   = fscanf(metadata_fid, '%d', 1);
    htree.minsize   = fscanf(metadata_fid, '%d', 1);
    htree.nnode     = fscanf(metadata_fid, '%d', 1);
    htree.root      = fscanf(metadata_fid, '%d', 1) + 1;
    htree.nlevel    = fscanf(metadata_fid, '%d', 1);
    h2mat.is_HSS    = fscanf(metadata_fid, '%d', 1);
    h2mat.minlvl    = fscanf(metadata_fid, '%d', 1) + 1;
    h2mat.n_near    = fscanf(metadata_fid, '%d', 1);
    h2mat.n_far     = fscanf(metadata_fid, '%d', 1);
    h2mat.par       = fscanf(metadata_fid, '%e', 1);
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
    for i = 1 : htree.nnode
        node = fscanf(metadata_fid, '%d', 1) + 1;
        lvl  = fscanf(metadata_fid, '%d', 1) + 1;
        htree.nodelvl(node)    = lvl;
        htree.cluster(node, 1) = fscanf(metadata_fid, '%d', 1) + 1;
        htree.cluster(node, 2) = fscanf(metadata_fid, '%d', 1) + 1;
        htree.level{lvl} = [htree.level{lvl} node];
        n_child = fscanf(metadata_fid, '%d', 1);
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

    %% 2. Binary data: input point coordinates and permutation array
    % The permutation array in H2Pack-Matlab is the inverse function of H2Pack permutation array
    coord0 = fread(binary_fid, [htree.dim, htree.n_point], 'double');
    coord0 = coord0';
    prem0  = fread(binary_fid, htree.n_point, 'int') + 1;
    perm   = zeros(htree.n_point, 1);
    coord  = coord0;
    for i = 1 : htree.n_point
        perm(prem0(i)) = i;
        coord(i, :) = coord0(prem0(i), :);
    end
    htree.permutation = perm;
    htree.coord = coord;

    %% 3. Metadata & binary data: U matrices
    h2mat.U = cell(htree.nnode, 1);
    for i = 1 : htree.nnode
        U_idx  = fscanf(metadata_fid, '%d', 1) + 1;
        U_nrow = fscanf(metadata_fid, '%d', 1);
        U_ncol = fscanf(metadata_fid, '%d', 1);
        tmpU   = fread(binary_fid, [U_ncol, U_nrow], 'double');
        h2mat.U{U_idx} = tmpU';
    end

    %% 4. Metadata & binary data: B matrices
    h2mat.far = zeros(h2mat.n_far, 2);
    h2mat.B   = cell(htree.nnode, htree.nnode);
    for i = 1 : h2mat.n_far
        node0  = fscanf(metadata_fid, '%d', 1) + 1;
        node1  = fscanf(metadata_fid, '%d', 1) + 1;
        B_nrow = fscanf(metadata_fid, '%d', 1);
        B_ncol = fscanf(metadata_fid, '%d', 1);
        tmpB   = fread(binary_fid, [B_ncol, B_nrow], 'double');
        h2mat.far(i, 1) = node0;
        h2mat.far(i, 2) = node1;
        h2mat.B{node0, node1} = tmpB';
    end

    %% 5. Metadata & binary data: D matrices
    h2mat.near = zeros(h2mat.n_near, 2);
    h2mat.D    = cell(htree.nnode, htree.nnode);
    for i = 1 : htree.n_leaf
        node   = fscanf(metadata_fid, '%d', 1) + 1;
        node   = fscanf(metadata_fid, '%d', 1) + 1;
        D_nrow = fscanf(metadata_fid, '%d', 1);
        D_ncol = fscanf(metadata_fid, '%d', 1);
        tmpD   = fread(binary_fid, [D_ncol, D_nrow], 'double');
        h2mat.D{node, node} = tmpD';
    end
    for i = 1 : h2mat.n_near
        node0  = fscanf(metadata_fid, '%d', 1) + 1;
        node1  = fscanf(metadata_fid, '%d', 1) + 1;
        D_nrow = fscanf(metadata_fid, '%d', 1);
        D_ncol = fscanf(metadata_fid, '%d', 1);
        tmpD   = fread(binary_fid, [D_ncol, D_nrow], 'double');
        h2mat.near(i, 1) = node0;
        h2mat.near(i, 2) = node1;
        h2mat.D{node0, node1} = tmpD';
    end

    %% 6. Metadata: H2Pack dedicated part: skeleton point indices
    has_skel = fscanf(metadata_fid, '%d', 1);
    h2mat.I = cell(htree.nnode, 1);
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
