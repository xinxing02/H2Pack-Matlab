function Yp = pp_numerical_selection_nlayer(kernel, L1, L2, L3, dim, tol)

    %   A special case for selecting proxy points only for H2 matrix
    %   with the two domains specified as 
    %       X = [-L1/2, L1/2]
    %       Y = [-L3/2, L3/2] \ [-L2/2, L2/2]

    %   kernel dimension
    kdim = size(kernel(randn(2, dim)), 1)/2;

    %   machine precision
    if isreal(kernel(randn(2, dim)))
       eps_machine = eps('double');
    else
       eps_machine = 100 * eps('double');
    end

    %   initial setup
    unit_box = [-ones(dim, 1), ones(dim, 1)];
    semi_L1 = 0.5*L1;
    semi_L2 = 0.5*L2;
    semi_L3 = 0.5*L3;
    box1  = semi_L1 * unit_box;
    domain1 = @(coord)domain_box(coord, box1);
    n_ring = 2*round((semi_L3 - semi_L2)/L1);
    %fprintf("%d annulus initalized\n", n_ring);
    
    %   Step 1: initial points in domain X
    Nx = 2000;    
    Ny = 5000;  % per ring
    X0 = sample_random(domain1, Nx, box1);
    
    %   Step 2: selection of Xp based on randomized compression of K(X0,
    %   Y0)
    %   Step 2.1: random multiplication without storing the kernel
    tmpA = zeros(kdim*Nx, kdim*Nx);
    for i = 1 : n_ring
        semi_Lmid1 = semi_L2 + (i-1)/n_ring*(semi_L3 - semi_L2);
        semi_Lmid2 = semi_L2 + i/n_ring    *(semi_L3 - semi_L2);
        hole_mid  = semi_Lmid1 * unit_box;
        box_mid   = semi_Lmid2 * unit_box;
        domain_mid = @(coord)domain_box_hole(coord, box_mid, hole_mid);
        Y0 = sample_random(domain_mid, Ny, box_mid);
        tmpA = tmpA + kernel({X0, Y0}) * randn(kdim*Ny, kdim*Nx);
    end
    
    %   Step 2.2: normlize columns and obtain Xp
    colnorm = sqrt(sum(tmpA.^2, 1));
    tmpA = bsxfun(@times, tmpA, 1./colnorm);     
    if (kdim == 1)
        [~, J] = ID(tmpA, 'reltol', tol);
    else
        [~, J] = IDdim(tmpA, kdim, 'reltol', tol);
        J = J(kdim:kdim:end)/kdim;
    end
    Xp = X0(J, :);

    %   Step 3: select the proxy points Yp in multiple steps. 
    Yp = zeros(0, dim);
    for i = 1 : n_ring
        semi_Lmid1 = semi_L2 + (i-1)/n_ring*(semi_L3 - semi_L2);
        semi_Lmid2 = semi_L2 + i/n_ring    *(semi_L3 - semi_L2);
        hole_mid  = semi_Lmid1 * unit_box;
        box_mid   = semi_Lmid2 * unit_box;
        domain_mid = @(coord)domain_box_hole(coord, box_mid, hole_mid);
        Y0 = sample_random(domain_mid, Ny, box_mid);
        tmpY = [Yp; Y0]; % the order of Yp and Y0 might matter
        A = kernel({Xp, tmpY});
        if kdim == 1
            [~, J] = ID(A', 'reltol', eps_machine);
        else
            [~, J] = IDdim(A', kdim, 'reltol', eps_machine);
            J = J(kdim:kdim:end)/kdim;
        end
        Yp = tmpY(J, :);
    end
    
    
%    fprintf("Initial selection gives %d proxy points\n", size(Yp, 1));  
%     %   Step 4: posterior correction for the selection Yp
%     %   Note by Xin: it seems that this correction doesn't have significant effect on the accuracy improments. 
%     Nx = 1000;    
%     Ny = 2000;
%     X0 = sample_random(domain1, Nx, box1);
%     fprintf("Posterior error check by spliting the far field into %d annulus\n", n_ring);
%     for i = 1 : n_ring
%         semi_Lmid1 = semi_L2 + (i-1)/n_ring*(semi_L3 - semi_L2);
%         semi_Lmid2 = semi_L2 + i/n_ring    *(semi_L3 - semi_L2);
%         hole_mid  = semi_Lmid1 * unit_box;
%         box_mid   = semi_Lmid2 * unit_box;
%         domain_mid = @(coord)domain_box_hole(coord, box_mid, hole_mid);
%         Y0 = sample_random(domain_mid, Ny, box_mid);
%         A = kernel({X0, Y0});
%         proxyA = kernel({X0, Yp});
%         [U, J] = ID(proxyA, 'reltol', tol);
%         normA = norm(A, 'fro');
%         if (normA < sqrt(Nx*Ny) * eps_machine) || (norm(A - U * A(J, :), 'fro') / normA < 2 * tol)
%             clear A proxyA Y0
%             fprintf("%dth annulus check: Yp works perfectly\n", i);
%         else
%             tmpY = [Y0; Yp];
%             A = kernel({Xp, tmpY});
%             [~, J] = ID(A', 'reltol', eps_machine);
%             Yp = tmpY(J, :);
%             clear A proxyA Y0
%             fprintf("%dth annulus check: large error detected and %d proxy points reselected\n", i,  size(Yp, 1));
%         end
%     end
end
