function logdet = HSS_ULV_logdet(ulvfactor, htree)


%   Use ULV decomposition of an HSS matrix to calculate the determinant

if isempty(ulvfactor)
    disp('ULV factors are not available');
    return
end


logdet = 0;

if ulvfactor.LU == true
    %   ULV decomposition based on LU:  A = L * U 
    Lc = ulvfactor.Lc;
    
    %   Splitted calculation   
    level = htree.level;
    nlevel = length(level);    
    
    for i = 1:nlevel
        for j = 1 : length(level{i})
            node = level{i}(j);      
            logdet = logdet + log(abs(det(Lc{node})));
        end
    end
    return ;
end


if ulvfactor.Chol == true
    L = ulvfactor.L;    
    level = htree.level;
    nlevel = length(level);    

    for i = 1:nlevel
        for j = 1 : length(level{i})
            node = level{i}(j);        
            tmp = sum(log(diag(L{node})));
            logdet = logdet + 2*tmp;
        end
    end    
end

end