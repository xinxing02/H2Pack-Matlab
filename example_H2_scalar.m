%%  Kernel Function/Matrix Definition
% kernel = @(coord)multiquadric(coord, 1/2, 100);
% kernel = @(coord)reciprocal(coord, 1);
kernel = @(coord)logfun(coord);
% kernel = @(coord)exponential(coord, 0.1);
% kernel = @(coord)gaussian(coord, 0.1);

%%  Point generation
npts = 80000;
dim = 2;
coord = npts^(1/3)*rand(npts, dim);
coord = bsxfun(@rdivide, coord, sum(coord.^2, 2).^(1/2));

%%  STEP1: Hierarchical partitioning 
minSize = 200;
htree = hierarchical_partition(coord, minSize, dim);
clear coord; % htree has a reorder copy of coord

%%  STEP2: Proxy Point Selection
alpha =  1;
reltol = 1e-5;
%   Scheme 1: proxy surface method, only works for kernel from
%   potential theory;
Yp = H2__ProxyPoint_Surface(htree, dim, alpha, 600);
%   Scheme 2: numerical selection, works for general kernels, using a
%   heuristic multi-step qr approach. 
% Yp = H2__ProxyPoint_QR_nlayer(kernel, htree, alpha, reltol);

%%  STEP3: H2 matrix construction using proxy points
JIT_flag = true; 
h2mat = Mat2H2_ID_Proxy(kernel, htree, Yp, 'reltol', reltol, alpha, JIT_flag);

%%  STEP4: H2 matrix-vector multiplications
x = randn(npts, 10);
tic
u_h2 = H2_matvec(h2mat, htree, x);
toc

%   error checking
idx = randperm(npts, 1000);
u_exact = kernel({htree.coord(idx, :), htree.coord}) * x;
err = sqrt(sum((u_h2(idx, :)-u_exact).^2, 1) ) ./ sqrt(sum(u_exact.^2, 1));

fprintf("min/mean/max relative errors for 10 matvecs:\n%.3e,%.3e,%.3e\n", ...
    min(err), mean(err), max(err));

%%  STEP5: H2 sub-matvec tests
nrow = 1000;
ncol = 112;
rowidx = randperm(npts, nrow);
colidx = randperm(npts, ncol);
x = randn(ncol, 10);
tic
u_h2 = H2_matvec_sub(h2mat, htree, x, rowidx, colidx);
toc

%   error checking
u_exact = kernel({htree.coord(rowidx, :), htree.coord(colidx, :)}) * x;
err = sqrt(sum((u_h2-u_exact).^2, 1) ) ./ sqrt(sum(u_exact.^2, 1));

fprintf("min/mean/max relative errors for 10 matvecs:\n%.3e,%.3e,%.3e\n", ...
    min(err), mean(err), max(err));


%%
%   Direct H2 matrix construction methods
%   WARNING: for large matrices, e.g., 40000*40000
%            it costs extremely large amount of storage and time
%
%$  METHOD 1: given a dense matrix
alpha = 1;
A = kernel(htree.coord);
h2mat = Mat2H2_ID(A, htree, 'reltol', 1e-6, alpha);

%%  METHOD 2: matrix-free, kernel blocks dynamically computed
alpha =  1;
JIT_flag = true;
h2mat = Mat2H2_ID_MatrixFree(kernel, htree, 'reltol', 1e-6, alpha, JIT_flag);

