function u = HSS_ULV_matvec(ulvfactor, htree, vec, flag)

%   Use ULV decomposition of an HSS matrix to multiply a vector

if isempty(ulvfactor)
    disp('ULV factors are not available');
    return
end

if ulvfactor.LU == true
    %   ULV decomposition based on LU:  A = L * U 
    %       if flag == 'L', u = L * vec;
    %       if flag == 'U', u = U * vec;
    %       if flag empty,  u = A * vec;
    
    if nargin < 4
        midu = HSS_ULV_matvec(ulvfactor, htree, vec, 'U');
        u    = HSS_ULV_matvec(ulvfactor, htree, midu, 'L');
        return 
    end
        
    Q = ulvfactor.Q;
    Lr = ulvfactor.Lr;
    Lc = ulvfactor.Lc;
    Idx = ulvfactor.Idx;
    
    level = htree.level;
    nlevel = length(level);    
    u = vec;
    
    if strcmp(flag, 'L')
        for i = 1 : nlevel
            for j = 1 : length(level{i})
                node = level{i}(j);                        
                u(Idx{node}, :) = Q{node} * (Lr{node} * u(Idx{node}, :));
            end
        end
    else
        for i = nlevel : -1 : 1
            for j = 1 : length(level{i})
                node = level{i}(j);                        
                u(Idx{node}, :) = Lc{node} * (Q{node}' * u(Idx{node}, :));
            end
        end
    end
    return ;
end

if ulvfactor.Chol == true
    %   ULV decomposition based on Chol:  A = L * L' 
    %       if flag == 'LT', u = L'* vec;
    %       if flag == 'L',  u = L * vec;
    %       if flag empty,   u = A * vec;
    
    if nargin < 4
        midu = HSS_ULV_matvec(ulvfactor, htree, vec, 'LT');
        u    = HSS_ULV_matvec(ulvfactor, htree, midu, 'L');
        return 
    end
    
    Q = ulvfactor.Q;
    L = ulvfactor.L;
    Idx = ulvfactor.Idx;
    
    level = htree.level;
    nlevel = length(level);    
    u = vec;
    
    if strcmp(flag, 'L')
        for i = 1 : nlevel
            for j = 1 : length(level{i})
                node = level{i}(j);                        
                u(Idx{node}, :) = Q{node} * (L{node} * u(Idx{node}, :));
            end
        end
    else
        for i = nlevel : -1 : 1
            for j = 1 : length(level{i})
                node = level{i}(j);                        
                u(Idx{node}, :) = L{node}' * (Q{node}' * u(Idx{node}, :));
            end
        end
    end
    return 
end