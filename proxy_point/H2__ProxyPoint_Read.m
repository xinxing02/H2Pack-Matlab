function Yp = H2__ProxyPoint_Read(filename, kernel, htree, tol)
%
%   Read the proxy points from the precomputed file. 
%   Note: if the corresponding box doesn't exist, compute the proxy points
%   and update the file as well. \
%
    

    
%     fprintf(fileID, dim, tol, alpha, minL, nlevel);

    
    %   partition basic info
    level  = htree.level;
    enbox  = htree.enbox;
    nlevel = length(level);
    Yp     = cell(nlevel, 1);
    L = enbox{htree.root}(2, :);
    dim = length(L);

    %   parse the stored info
    fileID = fopen(filename, 'r');
    para = sscanf(fgetl(fileID), '%d %e %d %f %d\n');
    numpp = sscanf(fgetl(fileID), '%d');
    linenum_pp = cumsum(numpp)'; 
    linenum_pp = [3, 3 + linenum_pp];  % the beginning line of the kth proxy point set
    dim0 = para(1);
    tol0 = para(2);
    nlayer0 = para(3);
    minL0 = para(4);
    nbox0 = para(5);
    if (dim0 ~= dim || abs(round(log2(L(1)/minL0)) - log2(L(1)/minL0)) > 1e-10 || tol0 > tol)
        disp 'An inconsistent proxy point file'
        Yp = {};
        return;
    end
    
    %   identify the corresponding boxes
    rootL = enbox{level{1}(1)}(2, 1);
    root_idx = round( log2(rootL/minL0) ) + 1;
    box_idx = root_idx - (2:(nlevel-1));  %    the boxes that we actually needed
                                        %    the file has boxes of idx 1,
                                        %    2, .., nbox
    
    %   select the proxy points
    for i = 1 : length(box_idx)
        id = box_idx(i);
        if (id >= 1) && (id <=nbox0)
            %   read lines linenum_pp(id) ~ linenum_pp(id+1)-1
            frewind(fileID);
            tmpYp = textscan(fileID, '%f', numpp(id)*dim, 'HeaderLines',linenum_pp(id)-1);
            Yp{i+2} = reshape(tmpYp{1}, dim, [])';
        else
            disp 'proxy points not found in the provided file';
            disp 'dynamically constructing the proxy points ...';
            L1 = enbox{level{i+2}(1)}(2, 1); 
            L2 = (1 + 2) * L1;
            L3 = min((1 + nlayer0)*L1, 2*rootL - L1);
            tmp_Yp = pp_numerical_selection_nlayer(kernel, L1, L2, L3, dim, tol);
            if isempty(tmp_Yp)
                Yp{i+2} = zeros(0, dim);
            else
                Yp{i+2} = tmp_Yp;
            end
        end
    end
    
    %   update the file if the smallest box < L0 or the largest box > L0*2^(nlvl0-1)
    if box_idx(end) < 1 || box_idx(1) > nbox0
        disp 'updating new proxy points into the same file'
        minL_new = minL0 * 2^(min(box_idx(end), 1) - 1);
        nbox_new = max(nbox0, box_idx(1)) - min(box_idx(end), 1) + 1;
        Ypp = cell(nbox_new, 1);
        numpp_new = zeros(nbox_new, 1);
        for i = 1 : nbox_new
            L1 = minL_new * 2^(i-1);
            old_idx = round(log2(L1/minL0)) + 1;
            if (old_idx >= 1) && (old_idx <=nbox0)
                %   read lines linenum_pp(id) ~ linenum_pp(id+1)-1
                frewind(fileID);
                tmpYp = textscan(fileID, '%f', numpp(old_idx)*dim, 'HeaderLines',linenum_pp(old_idx)-1);
                Ypp{i} = reshape(tmpYp{1}, dim, [])';
                numpp_new(i) = numpp(old_idx);
            else 
                L2 = (1 + 2) * L1;
                L3 = (1 + nlayer0)*L1;
                tmp_Yp = pp_numerical_selection_nlayer(kernel, L1, L2, L3, dim, tol0);
                if isempty(tmp_Yp)
                    Ypp{i} = zeros(0, dim);
                    numpp_new(i)=0;
                else
                    Ypp{i} = tmp_Yp;
                    numpp_new(i) = size(tmp_Yp, 1);
                end
            end
        end
        fclose(fileID);
        %%   store in file. 
        fileID = fopen(filename,'w');

        %   line 1: parameters
        fprintf(fileID, '%d %.3e %d %.12f %d\n', dim, tol0, nlayer0, minL_new, nbox_new);

        %   line 2: number of proxy points
        for i = 1 : nbox_new-1
            fprintf(fileID, '%d ', numpp_new(i));
        end
        fprintf(fileID, '%d\n', numpp_new(nbox_new));

        %   line 3-: proxy points
        for i = 1 : nbox_new
            Yp0 = Ypp{i};
            for j = 1 : numpp_new(i)
                for k = 1 : dim-1
                    fprintf(fileID, '%f ', Yp0(j, k));
                end
                fprintf(fileID, '%f\n', Yp0(j, end));
            end
        end
        fclose(fileID);
        
    end
        