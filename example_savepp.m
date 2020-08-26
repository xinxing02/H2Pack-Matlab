%%  Kernel Function, Relative Tol, Point dimension
kernel = @(coord)multiquadric(coord, 1/2, 100);
% kernel = @(coord)reciprocal(coord, 1);
% kernel = @(coord)exponential(coord, 0.1);
% kernel = @(coord)gaussian(coord, 0.1);
reltol = 1e-5;
dim = 3;

%%  Precompute and save proxy points
filename = './tmp_pp.dat';
minL = 1;   %   can initialize with better value if assuming a sample coord is available.
nbox = 6;
nlayer = 8;
precompute_proxy_point(filename, kernel, reltol, dim, minL, nbox, nlayer);

%%  Point generation
npts = 320000;
coord = npts^(1/3)*rand(npts, dim);

%%  STEP1: Hierarchical partitioning 
minSize = 200;
box = rootboxsize_proxy_point(coord, filename);
htree = hierarchical_partition(coord, minSize, dim, box);
clear coord; % htree has a reorder copy of coord

%%  STEP2: Proxy Point Selection
Yp = H2__ProxyPoint_Read(filename, kernel, htree, reltol);

%%  STEP3: H2 matrix construction using proxy points
alpha = 1;
JIT_flag = true; 
h2mat = Mat2H2_ID_Proxy(kernel, htree, Yp, 'reltol', reltol, alpha, JIT_flag);

%%  STEP4: H2 matrix-vector multiplications
x = randn(npts, 10);
u_h2 = H2_matvec(h2mat, htree, x);

%   error checking
idx = randperm(npts, 1000);
u_exact = kernel({htree.coord(idx, :), htree.coord}) * x;
err = sqrt(sum((u_h2(idx, :)-u_exact).^2, 1) ) ./ sqrt(sum(u_exact.^2, 1));

fprintf("min/mean/max relative errors for 10 matvecs:\n%.3e,%.3e,%.3e\n", ...
    min(err), mean(err), max(err));

