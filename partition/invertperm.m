function r = invertperm(p)
%
%   Invert the permutation vector
%   
    n = length(p);
    r = zeros(n,1);
    r(p) = (1:n)';
end