function box = rootboxsize_proxy_point(coord, filename)
%
%   Decide the box size according to the stored proxy point configuration
%   

fileID = fopen(filename, 'r');
% fprintf(fileID, '%d %.3e %d %.12f %d\n', dim, tol, nlayer, minL, nlevel);
para = sscanf(fgetl(fileID), '%d %e %d %f %d\n');
dim = para(1);
minL = para(4);

center = mean(coord, 1);
rela_coord = abs(bsxfun(@minus, coord, center));
L = max(rela_coord(:)); % half of the rootbox edge length
k = ceil(log2(L / minL));
L = minL * 2^k;         % half of the rootbox edge length
box(1, :) = center - L * ones(1, dim);  % corner of the box
box(2, :) = 2 * L * ones(1, dim);       % edgelength of the box.
fclose(fileID);

