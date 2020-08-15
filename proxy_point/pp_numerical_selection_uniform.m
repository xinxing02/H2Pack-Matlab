function Yp = pp_numerical_selection_uniform(kernel, domain1, box1, domain2, box2, tol, Nx, Ny)


if nargin <= 7 
    Nx = 2000;    
    Ny = 10000;
end

%   point dimension
pdim = size(box1, 1);
kdim = size(kernel(randn(2, pdim)),1)/2;

%   set up the machine precision
if isreal(kernel(randn(2, pdim)))
    eps_machine = eps('double');
else
    eps_machine = 100 * eps('double');
end

%   STEP 1: initial points in domain X and Y
X0 = sample_random(domain1, Nx, box1);
Y0 = sample_random(domain2, Ny, box2);
A = kernel({X0, Y0});

%   STEP 2: select Xp via randomized QR 
tmpA = A * randn(size(A,2), size(A,1));
colnorm = sqrt(sum(tmpA.^2, 1));
tmpA = bsxfun(@times, tmpA, 1./colnorm);     
if kdim == 1
    [~, J] = ID(tmpA, 'reltol', tol);
else
    [~, J] = IDdim(tmpA, kdim, 'reltol', tol);
end
Nx = min(Nx, length(J)/kdim);
% Xp = X0(J, :);

%   STEP 3: select proxy points Yp next. 
if kdim == 1
    [~, J] = ID(A(J, :)', 'reltol', eps_machine);
else
    [~, J] = IDdim(A(J, :)', kdim, 'reltol', eps_machine);
    J = J(kdim:kdim:end)/kdim;
end
Yp = Y0(J, :);

end
