function [err, blcnorm] = block__blockwise_error_func(kernel, coord, htree, blocklist, I, exU)
cluster = htree.cluster;
nodelvl = htree.nodelvl;
err = 0;
blcnorm = 0;
for i = 1 : size(blocklist, 1)
    c1 = blocklist(i,1);
    c2 = blocklist(i,2);
    idx1 = cluster(c1, 1):cluster(c1,2);
    idx2 = cluster(c2, 1):cluster(c2,2);
    
    %   Special case for extremely large sub-blocks
    if length(idx1) > 1e4
        [suberr, subblcnorm] = large_subblock(c1, c2, ceil(length(idx1)/5e3));
        err = err + suberr^2;
        blcnorm = blcnorm + subblcnorm^2;
    else        
        tmpA = kernel({coord(idx1,:), coord(idx2,:)}); 
        blcnorm = blcnorm + norm(tmpA, 'fro')^2;
        if isempty(I{c1}) || isempty(I{c2}) || isempty(exU{c1}) || isempty(exU{c2})
            disp 'given block not exist in the H2 representation';
            return ;
        end
        if nodelvl(c1) == nodelvl(c2)
            tmpAij = kernel({coord(I{c1},:), coord(I{c2},:)}); 
            err = err + norm(tmpA - exU{c1} * tmpAij * (exU{c2})', 'fro')^2;

    %         fprintf('block %d * %d has error ratio %e\n', c1, c2, ...
    %         norm(tmpA - exU{c1} * tmpAij * (exU{c2})', 'fro')/norm(tmpA, 'fro'));
        elseif nodelvl(c1) > nodelvl(c2)  
            % c2 is the leafnode at higher level. only compress on c1's side
            tmpAij = kernel({coord(I{c1},:), coord(idx2, :)});
            err = err + norm(tmpA - exU{c1} * tmpAij, 'fro')^2;

    %         fprintf('block %d * %d has error ratio %e\n', c1, c2, ...
    %         norm(tmpA - exU{c1} * tmpAij, 'fro')/norm(tmpA, 'fro'));
        else
            % c1 is the leafnode at higher level. only compress on c2's side
            tmpAij = kernel({coord(idx1, :), coord(I{c2}, :)});
            err = err + norm(tmpA - tmpAij * (exU{c2})', 'fro')^2;

    %         fprintf('block %d * %d has error ratio %e\n', c1, c2, ...
    %         norm(tmpA - tmpAij * (exU{c2})', 'fro')/norm(tmpA, 'fro'));
        end
    end
    clear tmpA tmpAij
end

blcnorm = sqrt(blcnorm);
err = sqrt(err);

function [suberr, subblcnorm] = large_subblock(c1, c2, nblk)
    %   basic check
    if isempty(I{c1}) || isempty(I{c2}) || isempty(exU{c1}) || isempty(exU{c2})
        disp 'given block not exist in the H2 representation';
        return ;
    end
    
    %   block cutting
    cut = round( (1:(nblk-1)) * (cluster(c1,2) - cluster(c1,1) + 1) / nblk );
    cut = [0, cut, cluster(c1,2)-cluster(c1,1)+1];
    subidxx = cell(nblk,1);
    for ii = 1 : nblk
        subidxx{ii} =  (cluster(c1,1)+cut(ii)):(cluster(c1,1)+cut(ii+1)-1);
    end
    cut = round( (1:(nblk-1)) * (cluster(c2,2) - cluster(c2,1) + 1) / nblk );
    cut = [0, cut, cluster(c2,2)-cluster(c2,1)+1];
    subidxy = cell(nblk,1);
    for ii = 1 : nblk
        subidxy{ii} =  (cluster(c2,1)+cut(ii)):(cluster(c2,1)+cut(ii+1)-1);
    end
    
    %   calculate error block-by-block
    if nodelvl(c1) == nodelvl(c2)
        tmpAij = kernel({coord(I{c1},:), coord(I{c2},:)}); 
        tmperr  = zeros(nblk*nblk,1);
        tmpnorm = zeros(nblk*nblk,1);
        parfor kk = 1 : (nblk*nblk)
            ii = mod(kk-1, nblk) + 1;
            jj = floor( (kk-1)/nblk ) + 1;
            tmpA = kernel({coord(subidxx{ii},:), coord(subidxy{jj}, :)});
            tmpnorm(kk) = norm(tmpA, 'fro')^2;        
            tmperr(kk)  = norm(tmpA - exU{c1}(subidxx{ii}-cluster(c1,1)+1, :) * ...
                    tmpAij * (exU{c2}(subidxy{jj}-cluster(c2,1)+1, :))', 'fro')^2; 
        end
        suberr = sqrt(sum(tmperr));
        subblcnorm = sqrt(sum(tmpnorm));
    elseif nodelvl(c1) > nodelvl(c2)  
        % c2 is the leafnode at higher level. only compress on c1's side
        tmpAij = kernel({coord(I{c1},:), coord(idx2,:)}); 
        tmperr  = zeros(nblk*nblk,1);
        tmpnorm = zeros(nblk*nblk,1);
        parfor kk = 1 : (nblk*nblk)
            ii = mod(kk-1, nblk) + 1;
            jj = floor( (kk-1)/nblk ) + 1;
            tmpA = kernel({coord(subidxx{ii},:), coord(subidxy{jj}, :)});
            tmpnorm(kk) = norm(tmpA, 'fro')^2;        
            tmperr(kk)  = norm(tmpA - exU{c1}(subidxx{ii}-cluster(c1,1)+1, :) * ...
                    tmpAij(:, subidxy{jj}-cluster(c2,1)+1), 'fro')^2; 
        end
        suberr = sqrt(sum(tmperr));
        subblcnorm = sqrt(sum(tmpnorm));
    else
        % c1 is the leafnode at higher level. only compress on c2's side
        tmpAij = kernel({coord(idx1,:), coord(I{c2},:)}); 
        tmperr  = zeros(nblk*nblk,1);
        tmpnorm = zeros(nblk*nblk,1);
        parfor kk = 1 : (nblk*nblk)
            ii = mod(kk-1, nblk) + 1;
            jj = floor( (kk-1)/nblk ) + 1;
            tmpA = kernel({coord(subidxx{ii},:), coord(subidxy{jj}, :)});
            tmpnorm(kk) = norm(tmpA, 'fro')^2;        
            tmperr(kk)  = norm(tmpA - tmpAij(subidxx{ii}-cluster(c1,1)+1, :) * ...
                (exU{c2}(subidxy{jj}-cluster(c2,1)+1, :))', 'fro')^2; 
        end
        suberr = sqrt(sum(tmperr));
        subblcnorm = sqrt(sum(tmpnorm));
    end
    
%     suberr = 0;
%     subblcnorm = 0;
%     tmpAij = kernel({coord(I{c1},:), coord(I{c2},:)}); 
%     for ii = 1 : nblk
%         for jj = 1 : nblk 
%             tmpA = kernel({coord(subidxx{ii},:), coord(subidxy{jj}, :)});
%             subblcnorm = subblcnorm + norm(tmpA, 'fro')^2;        
%             suberr = suberr + norm(tmpA - exU{c1}(subidxx{ii}-cluster(c1,1)+1, :) * ...
%                 tmpAij * (exU{c2}(subidxy{jj}-cluster(c2,1)+1, :))', 'fro')^2;     
%         end
%     end
%     suberr = sqrt(suberr);
%     subblcnorm = sqrt(subblcnorm);
end

end