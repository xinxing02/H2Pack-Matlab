kernel = @(coord)reciprocal(coord, 1);
[htree, h2mat] = read_H2_matrix_files('Coulomb_3D_1e-4.txt', 'Coulomb_3D_1e-4.bin');
h2mat.kernel = kernel;

%% Test AOT matvec
x       = randn(htree.n_point, 10);
u_h2    = H2_matvec(h2mat, htree, x);
idx     = randperm(htree.n_point, 1000);
u_exact = kernel({htree.coord(idx, :), htree.coord}) * x;
err     = sqrt(sum((u_h2(idx, :)-u_exact).^2, 1) ) ./ sqrt(sum(u_exact.^2, 1));

fprintf("min/mean/max relative errors for 10 AOT matvecs:\n%.3e,%.3e,%.3e\n", ...
    min(err), mean(err), max(err));

%% Test JIT matvec
h2mat.JIT = 1;
u_h2      = H2_matvec(h2mat, htree, x);
err       = sqrt(sum((u_h2(idx, :)-u_exact).^2, 1) ) ./ sqrt(sum(u_exact.^2, 1));
fprintf("min/mean/max relative errors for 10 JIT matvecs:\n%.3e,%.3e,%.3e\n", ...
    min(err), mean(err), max(err));
