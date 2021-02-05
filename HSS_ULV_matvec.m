function u = HSS_ULV_matvec(ulvfactor, htree, vec, kdim, flag)
%   Use ULV decomposition of an HSS matrix to multiply a vector

    %%  Permute the input vector 
    %vec(htree.permutation, :) = vec;
    vec0 = vec;
    npt = length(htree.permutation);
    for i = 1 : npt
        j     = htree.permutation(i);
        sidx0 = (i - 1) * kdim + 1;
        eidx0 =  i      * kdim;
        sidx1 = (j - 1) * kdim + 1;
        eidx1 =  j      * kdim;
        vec(sidx1 : eidx1, :) = vec0(sidx0 : eidx0, :);
    end

    if isempty(ulvfactor)
        disp('ULV factors are not available');
        return
    end

    if ulvfactor.LU == true
        if nargin < 5
            midu = HSS_ULV_matvec_LU(ulvfactor, htree, vec,  'U');
            u    = HSS_ULV_matvec_LU(ulvfactor, htree, midu, 'L');
        else
            u    = HSS_ULV_matvec_LU(ulvfactor, htree, vec,  flag);
        end
    end

    if ulvfactor.Chol == true
        if nargin < 5
            midu = HSS_ULV_matvec_Chol(ulvfactor, htree, vec,  'LT');
            u    = HSS_ULV_matvec_Chol(ulvfactor, htree, midu, 'L');
        else
            u    = HSS_ULV_matvec_Chol(ulvfactor, htree, vec,  flag);
        end
    end

    %%  Permute the output vector 
    %u = u(htree.permutation, :);
    u0 = u;
    npt = length(htree.permutation);
    for i = 1 : npt
        j     = htree.permutation(i);
        sidx0 = (i - 1) * kdim + 1;
        eidx0 =  i      * kdim;
        sidx1 = (j - 1) * kdim + 1;
        eidx1 =  j      * kdim;
        u(sidx0 : eidx0, :) = u0(sidx1 : eidx1, :);
    end
end

function u = HSS_ULV_matvec_LU(ulvfactor, htree, vec, flag)
%   ULV decomposition based on LU:  A = L * U 
%       if flag == 'L', u = L * vec;
%       if flag == 'U', u = U * vec;
%       if flag empty,  u = A * vec;

    Q      = ulvfactor.Q;
    Lr     = ulvfactor.Lr;
    Lc     = ulvfactor.Lc;
    Idx    = ulvfactor.Idx;
    level  = htree.level;
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
end

function u = HSS_ULV_matvec_Chol(ulvfactor, htree, vec, flag)
%   ULV decomposition based on Chol:  A = L * L' 
%       if flag == 'LT', u = L'* vec;
%       if flag == 'L',  u = L * vec;
%       if flag empty,   u = A * vec;

    Q      = ulvfactor.Q;
    L      = ulvfactor.L;
    Idx    = ulvfactor.Idx;
    level  = htree.level;
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
end