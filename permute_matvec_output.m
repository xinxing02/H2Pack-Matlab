function y = permute_matvec_output(perm, n_point, kdim, y0)
    ncol = size(y0, 2);
    y = zeros(size(y0));
    for i = 1 : ncol
        y0i = y0(:, i);
        yi  = zeros(n_point * kdim, 1);
        for j = 1 : n_point
            y0i_sidx = (perm(j) - 1) * kdim + 1;
            y0i_idx  = y0i_sidx : (perm(j) * kdim);
            yi_idx   = ((j - 1) * kdim + 1) : j * kdim;
            yi(yi_idx) = y0i(y0i_idx);
        end
        y(:, i) = yi;
    end
end