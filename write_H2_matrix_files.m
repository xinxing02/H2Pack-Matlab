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
    fprintf(metadata_fid, '%d\n', pt_dim);
    fprintf(metadata_fid, '%d\n', h2mat.kdim);
    fprintf(metadata_fid, '%d\n', n_point);
    fprintf(metadata_fid, '%d\n', h2mat.kdim * n_point);
    fprintf(metadata_fid, '%d\n', 1);  % If kernel matrix is symmetric
    fprintf(metadata_fid, '%d\n', htree.minsize);
    fprintf(metadata_fid, '%d\n', htree.nnode);
    fprintf(metadata_fid, '%d\n', htree.root - 1);
    fprintf(metadata_fid, '%d\n', htree.nlevel);
    fprintf(metadata_fid, '%d\n', is_HSS);
    fprintf(metadata_fid, '%d\n', h2mat.minlvl - 1);
    fprintf(metadata_fid, '%d\n', n_near);
    fprintf(metadata_fid, '%d\n', n_far);
    fprintf(metadata_fid, '%e\n', h2mat.par(1));
    for i = 1 : htree.nnode
        fprintf(metadata_fid, '%6d ', i - 1);
        fprintf(metadata_fid, '%2d ', htree.nodelvl(i) - 1);
        fprintf(metadata_fid, '%8d %8d ', htree.cluster(i, 1) - 1, htree.cluster(i, 2) - 1);
        node_i_childs = htree.children(i, :);
        n_child = sum(~isnan(node_i_childs));
        fprintf(metadata_fid, '%2d ', n_child);
        for j = 1 : n_child
            fprintf(metadata_fid, '%8d ', node_i_childs(j) - 1);
        end
        for j = n_child + 1 : max_child
            fprintf(metadata_fid, '-1 ');
        end
        fprintf(metadata_fid, '\n');
    end

    %% 2. Binary data: input point coordinates and permutation array
    % The permutation array in H2Pack-Matlab is the inverse function of H2Pack permutation array
    coord0 = zeros(n_point, pt_dim);
    perm   = htree.permutation;
    perm1  = zeros(size(perm));
    for i = 1 : n_point
        coord0(i, :) = htree.coord(perm(i), :);
        perm1(perm(i)) = i - 1;
    end
    coord0 = coord0';
    fwrite(binary_fid, coord0, 'double');
    fwrite(binary_fid, perm1, 'int');

    %% 3. Metadata & binary data: U matrices
    for i = 1 : htree.nnode
        tmpU = h2mat.U{i};
        [U_nrow, U_ncol] = size(tmpU);
        fprintf(metadata_fid, '%6d %5d %5d\n', i - 1, U_nrow, U_ncol);
        fwrite(binary_fid, tmpU', 'double');
    end

    %% 4. Metadata & binary data: B matrices
    for i = 1 : n_far
        node0 = h2mat.far(i, 1);
        node1 = h2mat.far(i, 2);
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
        fprintf(metadata_fid, '%6d %6d %5d %5d\n', node0 - 1, node1 - 1, B_nrow, B_ncol);
        fwrite(binary_fid, tmpB', 'double');
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
        fprintf(metadata_fid, '%6d %6d %5d %5d\n', node - 1, node - 1, D_nrow, D_ncol);
        fwrite(binary_fid, tmpD', 'double');
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
        fprintf(metadata_fid, '%6d %6d %5d %5d\n', node0 - 1, node1 - 1, D_nrow, D_ncol);
        fwrite(binary_fid, tmpD', 'double');
    end

    %% 6. Metadata: H2Pack dedicated part: skeleton point indices
    fprintf(metadata_fid, '%d\n', 1);  % Has H2Pack skeleton point indices
    for i = 1 : htree.nnode
        tmpI = h2mat.I{i};
        fprintf(metadata_fid, '%6d %6d ', i - 1, length(tmpI));
        for j = 1 : length(tmpI)
            fprintf(metadata_fid, '%d ', tmpI(j) - 1);
        end
        fprintf(metadata_fid, '\n');
    end

    fclose(metadata_fid);
    fclose(binary_fid);
end