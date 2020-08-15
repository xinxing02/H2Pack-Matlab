function [R, p] = sRRQR_ID(A, f, type, par)
%
%   Strong Rank Revealing QR (SPECIAL VERSION, ONLY FOR ID)
%       A P = [Q1, Q2] * [R11, R12; 
%                           0, R22]
%   where R11 and R12 satisfy that (inv(R11) * R12) has entries
%   bounded by a pre-specified constant which is not less than 1. 
%   
%   Input: 
%       A, matrix, target matrix that is appoximated.
%       f, scalar, constant that bounds the entries of calculated (inv(R11) * R12)
%    type, string, be either "rank" or "tol" and specify the way to decide
%                  the dimension of R11, i.e., the truncation. 
%     par, scalar, the parameter for "rank" or "tol" defined by the type. 
%
%   Output: 
%       R = [R11, R12];
%       p, permutation indices;
%   
%   Reference: 
%       Gu, Ming, and Stanley C. Eisenstat. "Efficient algorithms for 
%       computing a strong rank-revealing QR factorization." SIAM Journal 
%       on Scientific Computing 17.4 (1996): 848-869.
%
%   Note: 
%       This strong RRQR algorithm is specifically for ID and thus doesn't
%       provide the orthogonal matrix Q 


%   given a fixed rank 
if (strcmp(type, 'rank'))
    [R, p] = sRRQR_ID_rank(A, f, par);
    return ;
end

%   given a fixed error threshold
if (strcmp(type, 'tol'))
    [R, p] = sRRQR_ID_tol(A, f, par);
    return ;
end

%   report error otherwise
fprintf('No parameter type of %s !\n', type)        
R = [];
p = [];

end