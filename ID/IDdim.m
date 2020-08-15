function varargout = IDdim(A, dim, type, par)

%   Interpolative decomposition using the pivoted QR decomposition
%   for fat matrix A with batch size "dim";
%
%   Output: 
%       A(P^{-1}, :) = [eye; E] * A(J, :);
%   where P and J contains dim-length consecutive indices.

%%  Dimension of the matrix A
[m, n] = size(A);

%%  IDdim with fixed rank
if (strcmp(type, 'rank'))
    par = min([par, m, n]);
    [R, p, k] = QRdim_c(A', dim, type, par);
end

%%  IDdim with absolute tolerance 
if (strcmp(type, 'abstol'))
    [R, p, k] = QRdim_c(A', dim, 'tol', par);
end

%%  IDdim with relative tolerance 
if (strcmp(type, 'tol') || strcmp(type, 'reltol'))
    tol = par * max(sum(A.^2, 2).^(1/2)); 
    [R, p, k] = QRdim_c(A', dim, 'tol', tol);    
end

%%  Process of upper-triangular matrix
if (k == 0)   
    R = zeros(0, size(A, 1));
    p = 1 : size(A, 1);
else
    R = R(1:k, :);
end

%%  r, J, P
r = size(R, 1);                 %   obtained rank. 
J = p(1:r);                     %   skeleton rows
P = [];  P(p) = 1:length(p);    %   permutation indices

%%  Extra step (not necessary mostly) : sorting the index subset J and modify permutation indices P correspondingly
[J, I] = sort(J, 'ascend');
invI = zeros(length(I),1);
invI(I) = 1 : length(I);
invI = [invI', (length(I)+1):length(P)];
P = invI(P);  

%%   Output
if nargout == 3    
    %   [U, J] = ID(A, type, par, f)
    
    %   Special case : rank equals 0 (mainly from 'tol')
    if (r == 0)
        P = 1 : m;
        J = ones(0, 1);
        E = zeros(m, 0);
        varargout{1} = P;
        varargout{2} = J(:);
        varargout{3} = E;
        return ;
    end
    
    %   normal case
    E = linsolve(R(1:r, 1:r), R(1:r, (r+1):end), struct('UT', true))';  %(inv(R11)*R12)'
    E = E(:, I);    %   reordering due to the sorting of J. 
    varargout{1} = P;
    varargout{2} = J(:);
    varargout{3} = E;
    return ;
elseif nargout == 2
    %   [P, J, E] = ID(A, type, par, f);
    
    %   Special case : rank equals 0 (mainly from 'tol')
    if (r == 0)
        J = ones(0, 1);
        U = zeros(m, 0);
        varargout{1} = U;
        varargout{2} = J(:);
        return ;
    end
    
    %   normal case
    E = linsolve(R(1:r, 1:r), R(1:r, (r+1):end), struct('UT', true))';
    E = E(:, I);
    U = [eye(size(E,2)); E];
    U = U(P, :);
    varargout{1} = U;
    varargout{2} = J(:);
    return ;
end
