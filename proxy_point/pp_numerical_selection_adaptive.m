function Yp = pp_numerical_selection_adaptive(kernel, L1, L2, L3, dim, tol)

    %   A special case for selecting proxy points only for H2 matrix
    %   with the two domains specified as 
    %       X = [-L1/2, L1/2]
    %       Y = [-L3/2, L3/2] \ [-L2/2, L2/2]


    unit_box = [-ones(dim, 1), ones(dim, 1)];
    semi_L1 = 0.5*L1;
    semi_L2 = 0.5*L2;
    semi_L3 = 0.5*L3;
    box1  = semi_L1 * unit_box;
    hole2 = semi_L2 * unit_box;
    box2  = semi_L3 * unit_box;
    domain1 = @(coord)domain_box(coord, box1);
    domain2 = @(coord)domain_box_hole(coord, box2, hole2);
    
    %   Step 1: initial points in domain X and Y
    Nx = 2000;    
    Ny = 10000;
    X0 = sample_random(domain1, Nx, box1);
    Y0 = sample_random(domain2, Ny, box2);
    A  = kernel({X0, Y0});

    %   Step 2: select Xp such that K(Xp, y) can approximate 
    %   function K(x, y) with any x in domain Y

    %   2.1 randomized low-rank approximation of K(X0, Y0)
    tmpA = A * randn(size(A,2), size(A,1));
    colnorm = sqrt(sum(tmpA.^2, 1));
    tmpA = bsxfun(@times, tmpA, 1./colnorm);     

    %   2.2 truncated QR decomposition of the product
    [~, R, e] = qr(tmpA', 0);
    Nx = min([Nx, find(abs(diag(R)) > tol*abs(R(1,1)), 1, 'last')]);  
    Xp = X0(e(1:Nx), :);
    
    %   set up the machine precision
    if isreal(A)
       eps_machine = eps('double');
    else
       eps_machine = 100 * eps('double');
    end

    %   Step 3: select the proxy points Yp next. 
    [~, J] = ID(A(e(1:Nx), :)', 'reltol', eps_machine);
    Yp = Y0(J, :);
    clear A
    fprintf("initial selection based on %d sample points obtains %d proxy points\n", size(X0, 1), size(Yp, 1));

    %   Step 4: posterior improvement of the selection Yp
    Nx = 1000;    
    Ny = 5000;
    X0 = sample_random(domain1, Nx, box1);
    kk = max(4, 2*ceil((semi_L3 - semi_L2)/L1));
    fprintf("Posterior error check by spliting the far field into %d annulus\n", kk);
    for i = 1 : kk
        semi_Lmid1 = semi_L2 + (i-1)/kk*(semi_L3 - semi_L2);
        semi_Lmid2 = semi_L2 + i/kk    *(semi_L3 - semi_L2);
        hold_mid  = semi_Lmid1 * unit_box;
        box_mid   = semi_Lmid2 * unit_box;
        domain_mid = @(coord)domain_box_hole(coord, box_mid, hold_mid);
        Y0 = sample_random(domain_mid, Ny, box_mid);
        A = kernel({X0, Y0});
        proxyA = kernel({X0, Yp});
        [U, J] = ID(proxyA, 'reltol', tol);
        normA = norm(A, 'fro');
        if (normA < sqrt(Nx*Ny) * eps_machine) || (norm(A - U * A(J, :), 'fro') / normA < 2 * tol)
            clear A proxyA Y0
            fprintf("%dth annulus check: Yp works perfectly\n", i);
        else
            tmpY = [Y0; Yp];
            A = kernel({Xp, tmpY});
            [~, J] = ID(A', 'reltol', eps_machine);
            Yp = tmpY(J, :);
            clear A
            fprintf("%dth annulus check: large error detected and %d proxy points reselected\n", i,  size(Yp, 1));
        end
    end
end
