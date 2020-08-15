function varargout = QRdim(A, dim, type, par)
%
%   block version of the QR decomposition
%

[m, n] = size(A);
mn = min([m,n]);
if strcmp(type, 'rank')
    par = min([par, m, n]);
end
colnorm = sum(A.^2, 1).^(1/2);
H = cell(mn, 1);
p = 1:n;
r = -1; 


for i = 1 : (mn/dim)        
    workidx = (i-1)*dim + (1:dim);

    %   calculate the column block's norm and find the pivot
    blknorm = arrayfun(@(j)norm(colnorm((j-1)*dim + (1:dim))), 1 : (n/dim));
    [val, ii] = max(blknorm);
%     fprintf('%d, %e\n', i, val);
    swapidx = (ii-1)*dim + (1:dim);

    %   check the termination
    if (strcmp(type, 'rank') && (i-1)*dim >= par)
        r = (i-1) * dim;
        break;
    end
    if (strcmp(type, 'tol') && val < dim^(1/2) * par)
        r = (i-1) * dim;
        break;
    end

    %   swap the data first
    colnorm([workidx, swapidx]) = colnorm([swapidx, workidx]);    
    p([workidx, swapidx]) = p([swapidx, workidx]);
    A(:, [workidx, swapidx]) = A(:, [swapidx, workidx]);

    %   orthogonalization 
    for j = workidx
        if (j == mn)            
            r = mn;
            break;
        end
        v = A(j:end, j);
        sgn = sign(v(1));
        if (sgn == 0); sgn = 1; end
        A(j, j) = - sgn * norm(v);
        A(j+1:end, j) = 0;
        v(1) = v(1) + sgn * norm(v);
        v = v / norm(v);
        H{j} = v;        
%         A(j:end, j+1:end) = A(j:end, j+1:end) - 2 * v * (v' * A(j:end, j+1:end));
        for k = (j+1) : n
            A(j:end, k) = A(j:end, k) - 2 * v * (v' * A(j:end, k));
        end                
        colnorm(i*dim+1:end) = (colnorm(i*dim+1:end).^2 - A(j, i*dim+1:end).^2).^(1/2);
    end        
    colnorm(workidx) = 0;
end



%   complete QR decomposition
if nargout == 3
    Q = eye(m);    
    for k = min([r,mn-1]) : -1 : 1
        Q(k:end, k:end) = Q(k:end, k:end) - 2 * H{k} * (H{k}' * Q(k:end, k:end));
    end
    Q = Q(:, 1:r);
    A = A(1:r, :);
    varargout{1} = Q;
    varargout{2} = A;
    varargout{3} = p;
else
    A = A(1:r, :);
    varargout{1} = A;
    varargout{2} = p;
end
     
% 
% 
% if mod(n, dim) > 0
%     Q = []; R = []; p = [];
%     return ;
% end
%    
% for i = 1 : (n/dim)
%     workidx = (i-1)*dim + (1:dim);
% 
%       calculate the column block's norm and find the pivot
%     blknorm = arrayfun(@(j)norm(colnorm((j-1)*dim + (1:dim))), 1 : (n/dim));
%     [val, ii] = max(blknorm);
%     fprintf('%d, %e\n', i, val);
%     swapidx = (ii-1)*dim + (1:dim);
% 
%       check the termination
%     if (strcmp(type, 'rank') && (i-1)*dim >= par)
%         r = (i-1) * dim;
%         break;
%     end
%     if (strcmp(type, 'tol') && val < dim^(1/2) * par)
%         r = (i-1) * dim;
%         break
%     end
%     
%       swap the data first
%     colnorm([workidx, swapidx]) = colnorm([swapidx, workidx]);
%     colnorm(workidx) = 0;
%     p([workidx, swapidx]) = p([swapidx, workidx]);
%     A(:, [workidx, swapidx]) = A(:, [swapidx, workidx]);
%     R(:, [workidx, swapidx]) = R(:, [swapidx, workidx]);
% 
%       orthogonalization 
%     [tmpQ, tmpR] = qr(A(:, workidx), 0);
%     R(workidx, workidx) = tmpR;    
%     R(workidx, i*dim+1:end) = tmpQ' * A(:, i*dim+1:end);
%     parfor k = (i*dim+1) : n
%         tmp = tmpQ' * A(:, k);
%         A(:, k)  = A(:, k) - tmpQ * tmp;
%     end
%     for k = (i*dim+1) : n
%         A(:, k) = A(:, k) - tmpQ * R(workidx, k);
%     end
%     A(:, i*dim+1:end) = A(:, i*dim+1:end) - tmpQ * R(workidx, i*dim+1:end);
%     Q(:, workidx) = tmpQ;
%     colnorm(i*dim+1:end) = (colnorm(i*dim+1:end).^2 - sum(R(workidx, i*dim+1:end).^2,1)).^(1/2);
% end

end
