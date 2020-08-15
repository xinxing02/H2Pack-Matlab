%%  Kernel Function/Matrix Definition
radii = 1;
eta   = 1/6/pi/radii;
kernel  = @(coord)rpy_mex(coord, radii, eta);
kdim = 3;

%%  Random generation of points
npts = 40000;
dim = 3;
coord = 8*npts^(1/dim) * [rand(npts, 2), zeros(npts, 1)];

%%  STEP1: Hierarchical partitioning 
minSize = 200;
htree = hierarchical_partition(coord, minSize, dim);
clear coord; % htree has a reorder copy of coord

%%  STEP2: Proxy Point Selection for far field interaction
alpha_pp =  1;
reltol = 1e-3;
%   Scheme 1: proxy surface method, only works for kernel from
%   potential theory;
% Yp = H2__ProxyPoint_Surface(htree, dim, alpha_pp, 600);
%   Scheme 2: numerical selection, works for general kernels, using a
%   heuristic multi-step qr approach. 
Yp = H2__ProxyPoint_QR_nlayer(kernel, htree, alpha_pp, reltol);

%%  STEP3: HSS matrix construction using Yp + nearfield blocks
JIT_flag = true;
hssmat = Mat2HSS_Hybrid(kernel, htree, Yp, 'reltol', reltol, alpha_pp, JIT_flag);

%%  STEP4: HSS matrix-vector multiplications
x = randn(kdim*npts, 10);
u_hss = H2_matvec(hssmat, htree, x);
%   error checking
idx = 1:600;
u_exact = kernel({htree.coord(idx, :), htree.coord}) * x;
err = sqrt(sum((u_hss(1:(600*kdim), :)-u_exact).^2, 1) ) ./ sqrt(sum(u_exact.^2, 1));
fprintf("min/mean/max relative errors in 10 matvec:\n%.3e,%.3e,%.3e\n", ...
    min(err), mean(err), max(err));

%%  STEP5: ULV decomposition 
%   NOTE: When the matrix itself is extremely illconditioned, ULV decomposition
%         and ULV solve can be unstable. Usually this system requires a diagonal 
%         shift.

% ulvfactor = HSS_ULV_Chol(hssmat, htree);  %for SPD HSS matrix 
ulvfactor = HSS_ULV_LU(hssmat, htree);  %for general symmetric HSS matrix

%%  STEP5.1: ULV solve
x = randn(kdim*npts, 10);
b = H2_matvec(hssmat, htree, x);
x_sol = HSS_ULV_solve(ulvfactor, htree, b);
err = sqrt(sum((x-x_sol).^2, 1) ) ./ sqrt(sum(x_sol.^2, 1));
fprintf("min/mean/max relative errors in 10 ULV-solves:\n%.3e,%.3e,%.3e\n", ...
    min(err), mean(err), max(err));

%%  STEP5.2: ULV matvec
x = randn(kdim*npts, 10);
b_hss = H2_matvec(hssmat, htree, x);
b_ulv = HSS_ULV_matvec(ulvfactor, htree, x);
err = sqrt(sum((b_hss - b_ulv).^2, 1) ) ./ sqrt(sum(b_hss.^2, 1));
fprintf("min/mean/max relative errors in 10 ULV-matvecs:\n%.3e,%.3e,%.3e\n", ...
    min(err), mean(err), max(err));


%%
%   Direct HSS matrix construction methods
%   WARNING: for large matrices, e.g., 40000*40000
%            it costs extremely large amount of storage and time
%%  METHOD 1: Direct HSS construction
alpha = 0;
A = kernel(htree.coord);
hssmat = Mat2H2_ID(A, htree, 'reltol', 1e-6, alpha);

%%  METHOD 2: Direct H2 construction (matrix free)
alpha =  0;
JIT_flag = true;
hssmat = Mat2H2_ID_MatrixFree(kernel, htree, 'reltol', 1e-6, alpha, JIT_flag);
