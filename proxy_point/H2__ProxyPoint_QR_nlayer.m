function Yp = H2__ProxyPoint_QR_nlayer(kernel, htree, alpha, tol)

    %   basic info
    level  = htree.level;
    enbox  = htree.enbox;
    nlevel = length(level);
    Yp     = cell(nlevel, 1);
    L = enbox{htree.root}(2, :);
    dim = length(L);

    %   selection of proxy points at each level
    maxL = enbox{level{1}(1)}(2, 1);
    for i = 3 : nlevel 
        %   always assume the enclosing boxes are cubic.
        L1 = enbox{level{i}(1)}(2, 1); 
        L2 = (1 + 2 * alpha) * L1;
        L3 = min((1 + 8 * alpha)*L1, 2*maxL - L1);
        tmp_Yp = pp_numerical_selection_nlayer(kernel, L1, L2, L3, dim, tol);
        if isempty(tmp_Yp)
            Yp{i} = zeros(0, dim);
        else
            Yp{i} = tmp_Yp;
        end
    end
