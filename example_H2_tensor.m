%%  Kernel Function/Matrix Definition
%   RPY kernel function is a 3*3 tensor kernel function from potential
%   theory and thus the proxy surface method works as well.
radii = 1;
eta   = 1/6/pi/radii;
kernel  = @(coord)rpy_mex(coord, radii, eta);
kdim = 3;

%%  Point Generation
npts = 30000;
dim = 3;
coord = 8*npts^(1/dim) * rand(npts, dim);

%%  Step1: Hierarchical partitioning 
minSize = 300;
htree = hierarchical_partition(coord, minSize, dim);

%%  STEP2: Proxy Point Selection
alpha =  1;
reltol = 1e-3;
%   Scheme 1: proxy surface method, only works for kernel from
%   potential theory;
% Yp = H2__ProxyPoint_Surface(htree, dim, alpha, 600);
%   Scheme 2: numerical selection, works for general kernels, using a
%   heuristic multi-step qr approach. 
Yp = H2__ProxyPoint_QR_nlayer(kernel, htree, alpha, reltol);

%%  STEP3: H2 matrix construction using proxy points
JIT_flag = true; 
h2mat = Mat2H2_ID_Proxy(kernel, htree, Yp, 'reltol', reltol, alpha, JIT_flag);

%%  STEP4: H2 matrix-vector multiplications
x = randn(kdim*npts, 10);
u_h2 = H2_matvec(h2mat, htree, x);
%   error checking at the first 1800 entries
idx = 1 : 600;
u_exact = kernel({htree.coord(idx, :), htree.coord}) * x;
err = sqrt(sum((u_h2(1:600*kdim, :)-u_exact).^2, 1) ) ./ sqrt(sum(u_exact.^2, 1));
fprintf("min/mean/max relative errors in 10 matvec:\n%.3e,%.3e,%.3e\n", ...
    min(err), mean(err), max(err));


%%
%   Direct H2 matrix construction methods
%   WARNING: for large matrices, e.g., 40000*40000
%            it costs extremely large amount of storage and time
%
%%  METHOD 1: given a dense matrix
%  WARNING: for large matrices, e.g., 40000*40000
%           it costs extremely large amount of storage and time
alpha = 1;
A = kernel(htree.coord);
h2mat = Mat2H2_ID(A, htree, 'reltol', 1e-6, alpha);

%%  METHOD 2: matrix-free, kernel blocks dynamically computed
alpha =  1;
JIT_flag = true;
h2mat = Mat2H2_ID_MatrixFree(kernel, htree, 'reltol', 1e-6, alpha, JIT_flag);