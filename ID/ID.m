function varargout = ID(A, type, par, f)

%   Interpolative Decomposition using Strong RRQR over rows of a target
%   matrix
%
%   Given an m*n matrix A, an rank-k ID approximation of A is of form
%       A = U * A(J, :)
%   where J is a row index subset of size k, and U is a m*k matrix with
%   entries bounded by the prespecified constant 'f'. A(J,:) and U are
%   usually called the skeleton and projection matrix, respectively. 
%   
%   Ref: 
%       Cheng, Hongwei, et al. "On the compression of low rank matrices." 
%       SIAM Journal on Scientific Computing 26.4 (2005): 1389-1404.
%
%   Input: 
%          A, target matrix
%       type, 'rank' or 'tol', the way to decide the rank of ID
%        par, parameter for the given 'type'
%          f, constant that bounds the entries of projectio matrix U in ID.
%   
%   Output: 
%           [U, J] = ID(A, type, par, f)
%                  ==> A \approx U * A(J, :)
%
%        [P, J, E] = ID(A, type, par, f);
%                  ==> P is a permuted row indices of A and E is a small
%                  matrix. They are used to construct U as
%                       tmp = [eye; E];
%                       U = tmp(P, :);
%                  and A \approx U * A(J, :)
%   Note: storing P and E is cheaper than storing U
%                 
%   Algorithm description: 
%       The ID of A depends on the strong RRQR of A' as follows. 
%           A' P = [Q1, Q2] [R11, R12; 
%                            0, R22]
%                \approx Q1 [R11, R12]
%                =  [A1', A1'*inv(R11)*R12];
%                = A1' * [eye, inv(R11)*R12];
%                = A1' * [eye, E'] (denoted as E')
%   Therefore, 
%           P*A = [eye; E] * A1; 
%   and denote U = P*[eye;E], A1 = A(J, :), then 
%             A \approx U * A(J, :);
%   Furthermore, with strong RRQR, entries of U are bounded by the given
%   parameter. 
%
%   The dimension of R11, i.e., length(J), is decided by the strong RRQR of
%   A' with the given type and par. 

%%  Initialization 
if nargin < 4 
    f = 2;
end
if f < 1
    f = 2;
end
[m, n] = size(A);

if n == 0
    r = 0;
else    
    %%  ID with fixed rank. 
    if (strcmp(type, 'rank'))
        %   modify rank if necessary.
        par = min([par, m, n]);
        %   construction
        [R, p] = sRRQR_ID(A', f, type, par);
        r = size(R, 1);                 %   obtained rank. 
        J = p(1:r);                     %   skeleton rows
        P = [];  P(p) = 1:length(p);    %   permutation indices
    end

    %%  ID with absolute tolerance
    if (strcmp(type, 'tol') || strcmp(type, 'abstol'))   
        [R, p] = sRRQR_ID(A', f, 'tol', par);
        r = size(R, 1);                 %   obtained rank. 
        J = p(1:r);                     %   skeleton rows
        P = [];  P(p) = 1:length(p);    %   permutation indices
    end


    %%  ID with relative tolerance
    if (strcmp(type, 'reltol'))
        %   calculate the abs tolerance
        tol = par * max(sum(A.^2, 2).^(1/2));
        %   construction
        [R, p] = sRRQR_ID(A', f, 'tol', tol);  
        r = size(R, 1);                 %   obtained rank. 
        J = p(1:r);                     %   skeleton rows
        P = [];  P(p) = 1:length(p);    %   permutation indices
    end

    %%   Extra step (not necessary mostly) : sorting the index subset J and modify permutation indices P correspondingly
    [J, I] = sort(J, 'ascend');
    invI = zeros(length(I),1);
    invI(I) = 1 : length(I);
    invI = [invI', (length(I)+1):length(P)];
    P = invI(P);  
end


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
