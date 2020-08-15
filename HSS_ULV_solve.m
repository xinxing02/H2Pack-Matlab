function u = HSS_ULV_solve(ulvfactor, htree, rhs, flag)

%   Use ULV decomposition of an HSS matrix to solve a linear system

if isempty(ulvfactor)
    disp('ULV factors are not available');
    return
end


if ulvfactor.LU == true
    %   ULV decomposition based on LU:  A = L * U 
    %       if flag == 'L', solve L * u = rhs;
    %       if flag == 'U', solve U * u = rhs;
    %       if flag empty,  solve A * u = rhs;
    
    if nargin < 4
        midu = HSS_ULV_solve(ulvfactor, htree, rhs, 'U');
        u    = HSS_ULV_solve(ulvfactor, htree, midu, 'L');
        return ;
    end
    
    Q = ulvfactor.Q;
    Lr = ulvfactor.Lr;
    Lc = ulvfactor.Lc;
    Idx = ulvfactor.Idx;
    
    %   Splitted calculation   
    level = htree.level;
    nlevel = length(level);    
    u = rhs;
    
    if strcmp(flag, 'L')
        for i = 1:nlevel
            for j = 1 : length(level{i})
                node = level{i}(j);                        
                u(Idx{node}, :) = Q{node} * (Lc{node}\u(Idx{node}, :));
            end
        end
    else
        for i = nlevel : -1 : 1
            for j = 1 : length(level{i})
                node = level{i}(j);                        
                u(Idx{node}, :) = Lr{node} \ (Q{node}'*u(Idx{node}, :));
            end
        end
    end
    return ;
end


if ulvfactor.Chol == true
    %   ULV decomposition based on Chol:  A = L * L' 
    %       if flag == 'LT', solve L'* u = rhs;
    %       if flag == 'L',  solve L * u = rhs;
    %       if flag empty,   solve A * u = rhs;
    
    if nargin < 4 
        midu = HSS_ULV_solve(ulvfactor, htree, rhs, 'LT');
        u    = HSS_ULV_solve(ulvfactor, htree, midu, 'L');
        return ;
    end
        
    Q = ulvfactor.Q;
    L = ulvfactor.L;
    r = ulvfactor.r;
    Idx = ulvfactor.Idx;
    
    %   Splitted calculation  
    level = htree.level;
    nlevel = length(level);    
    u = rhs;
    
    if strcmp(flag, 'L')
        for i = 1:nlevel
            for j = 1 : length(level{i})
                node = level{i}(j); 
                rr = r(node);   
                %   compute tmp = L{node}'\u(Idx{node}, :);
                tmp = u(Idx{node}, :);                
                tmp2 = tmp(rr+1:end, :) - L{node}(1:rr, rr+1:end)' * tmp(1:rr, :);
                tmp(rr+1:end, :) = linsolve(L{node}(rr+1:end, rr+1:end)', tmp2, struct('UT', true));
                u(Idx{node}, :) = Q{node} * tmp;
            end
        end
    else
        for i = nlevel : -1 : 1
            for j = 1 : length(level{i})
                node = level{i}(j);        
                rr = r(node);   
                tmp = Q{node}'*u(Idx{node}, :);
                %   compute tmp = L{node} \ tmp;
                x2 = linsolve(L{node}(rr+1:end, rr+1:end), tmp(rr+1:end, :), struct('LT', true));                
                u(Idx{node}, :) = [tmp(1:rr, :) - L{node}(1:rr, rr+1:end) * x2; x2];
            end
        end
    end
    return ;
end