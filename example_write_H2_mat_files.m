%%  Kernel Function/Matrix Definition
kernel = @(coord)reciprocal(coord, 1);

%%  Point generation
npts  = 40000;
dim   = 3;
coord = npts^(1/dim)*rand(npts, dim);

%%  STEP1: Hierarchical partitioning 
minsize = 400;
htree   = hierarchical_partition(coord, minsize, dim);
clear coord; % htree has a reorder copy of coord

test_H2  = 1;
test_HSS = (1 - test_H2);
if (test_H2 == 1)
    %%  STEP2: Proxy Point Selection
    alpha  = 1;
    reltol = 1e-4;
    %   Scheme 1: proxy surface method, only works for kernel from
    %   potential theory;
    % Yp = H2__ProxyPoint_Surface(htree, dim, alpha, 600);
    %   Scheme 2: numerical selection, works for general kernels, using a
    %   heuristic multi-step qr approach. 
    Yp = H2__ProxyPoint_QR_nlayer(kernel, htree, alpha, reltol);

    %%  STEP3: H2 matrix construction using proxy points
    JIT_flag = true; 
    h2mat = Mat2H2_ID_Proxy(kernel, htree, Yp, 'reltol', reltol, alpha, JIT_flag);

    %%  STEP4: Write H2 matrix data to file
    write_H2_matrix_files(htree, h2mat, 'Coulomb_3D_1e-4.txt', 'Coulomb_3D_1e-4.bin');
end
if (test_HSS == 1)
    %%  STEP2: Proxy Point Selection for far field interaction
    alpha_pp = 1;
    reltol   = 1e-4;
    %   Scheme 1: proxy surface method, only works for kernel from
    %   potential theory;
    % Yp = H2__ProxyPoint_Surface(htree, dim, alpha_pp, 600);
    %   Scheme 2: numerical selection, works for general kernels, using a
    %   heuristic multi-step qr approach. 
    Yp = H2__ProxyPoint_QR_nlayer(kernel, htree, alpha_pp, 1e-2 * reltol);

    %%  STEP3: HSS matrix construction using Yp + nearfield blocks
    JIT_flag = true;
    hssmat = Mat2HSS_Hybrid(kernel, htree, Yp, 'reltol', 1e-2 * reltol, alpha_pp, JIT_flag);

    %%  STEP4: Write H2 matrix data to file
    write_H2_matrix_files(htree, hssmat, 'Coulomb_3D_1e-4.txt', 'Coulomb_3D_1e-4.bin');
end