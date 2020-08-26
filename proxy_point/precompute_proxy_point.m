function precompute_proxy_point(filename, kernel, tol, dim, minL, nbox, nlayer)

%   
%   data structure of the file
%   line 1. dim,   point dimension
%           tol,   reltol
%           nlayer, number of layers in the annulus (not precisely)
%           minL,    minimum box edge length: 
%           nlvl,  number of boxes concerned
%   line 2. npp,   array of nlvl numbers, the number of proxy points for each
%   source box.
%   line 3-(npp(1)+2): proxy point coordinates for box 1
%   line (npp(1)+3)-(npp(1)+npp(2)+2): proxy point coordinates for box 2


%%   proxy point selection
Yp = cell(nbox, 1);
numpp = zeros(nbox, 1);
for i = 1 : nbox
    %   compute the proxy points for the box of size minL * 2^(i-1)
    L1 = minL * 2^(i-1);
    L2 = (1 + 2)*L1;
    L3 = (1 + nlayer) * L1;
    tmp_Yp = pp_numerical_selection_nlayer(kernel, L1, L2, L3, dim, tol);
    if isempty(tmp_Yp)
        Yp{i} = zeros(0, dim);
        numpp(i) = 0;
    else
        Yp{i} = tmp_Yp;
        numpp(i) = size(tmp_Yp, 1);
    end
end

%%   store in file. 
fileID = fopen(filename,'w');

%   line 1: parameters
fprintf(fileID, '%d %.3e %d %.12f %d\n', dim, tol, nlayer, minL, nbox);

%   line 2: number of proxy points
for i = 1 : nbox-1
    fprintf(fileID, '%d ', numpp(i));
end
fprintf(fileID, '%d\n', numpp(nbox));

%   line 3-: proxy points
for i = 1 : nbox
    Yp0 = Yp{i};
    for j = 1 : numpp(i)
        for k = 1 : dim-1
            fprintf(fileID, '%f ', Yp0(j, k));
        end
        fprintf(fileID, '%f\n', Yp0(j, end));
    end
end
fclose(fileID);
end
