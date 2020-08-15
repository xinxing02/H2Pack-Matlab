function [err, blcnorm] = block__blockwise_error(A, blocklist, mcluster, I, exU)
err = 0;
blcnorm = 0;
for i = 1 : length(blocklist)
    c1 = blocklist(i,1);
    c2 = blocklist(i,2);
    idx1 = mcluster(c1, 1):mcluster(c1,2);
    idx2 = mcluster(c2, 1):mcluster(c2,2);
    blcnorm = blcnorm + norm(A(idx1, idx2), 'fro')^2;
    if isempty(I{c1}) || isempty(I{c2}) || isempty(exU{c1}) || isempty(exU{c2})
        disp 'given block not exist in the H2 representation';
        return ;
    end
    err = err + norm(A(idx1, idx2) - exU{c1} * A(I{c1}, I{c2}) * (exU{c2})', 'fro')^2;
end

blcnorm = sqrt(blcnorm);
err = sqrt(err);

