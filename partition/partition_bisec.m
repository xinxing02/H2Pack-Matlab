function [parent, children, cluster, enbox, permutation, n] = ...
                        partition_bisec(coord, minSize, dim, box)
%   Hierarchical partitioning of the given points 'coord' in dim-D space.
%   Post-order tree is used
%
%   Input: 
%       coord,      N*3 matrix
%       minSize,    minimum size for the leaf nodes
%       dim,        dimension of the points. 
%       box,        2*dim matrix representing the box that encloses the
%                   point set 'coord'. 
%                   box(1, :) gives the coordinates of the corner of the box 
%                   with minimal x/y/z. 
%                   box(2, :) gives the edge length of the box along each
%                   dimension. 
%   
%   Output:
%       n:          number of nodes in the tree
%       parent:     n*1 vector, indicating the parent of the node. 
%       children:   n*2^dim matrix, indicating the children of the node.
%       cluster:    n*2 matrix, defining the indices at each cluster,
%                       [start-point, end-point]
%       enbox:      n*1 cell, each contains the box data that encloses the
%                       the associated cluster
%       permutation:n*1 vector, the permutation of the coordinates that makes 
%                       the clusters in order 

    N = size(coord,1);
    nchildren = power(2, dim);
    
    %%   Check dimension
    if size(coord,2) < dim
        disp 'Inconsistent coordinate dimension with the given dim'
        parent = []; children = []; cluster = []; enbox = {}; permutation = []; n = 0;
        return ;
    else
        coord(:,(dim+1):end) = [];
    end
    
    %%   Parse the minSize constraint
    if length(minSize) > 1
        minBox = minSize(2);
        minNum = minSize(1);
    else
        minBox = 0;
        minNum = minSize;
    end
    
    %%   Check the enclosing box
    if nargin < 4
        box = zeros(2, dim);
        center = mean(coord, 1);
        rela_coord = abs(bsxfun(@minus, coord, center));
        L = max(rela_coord(:)) + 1e-10;
        box(1, :) = center - L * ones(1, dim);
        box(2, :) = 2 * L * ones(1, dim);
        
%         box(1, :) = min(coord, [], 1);
%         box(2, :) = max(bsxfun(@minus, coord, box(1, :)), [], 1);
        if min(abs(box(2,:))) < eps('double')
            disp 'point lying in a lower-dimension'
            parent = []; children = []; cluster = []; permutation = []; n = 0;
            return ;
        else
            %   minor modification to avoid rounding error
            box(1,:) = box(1,:) - 200*eps('double');
            box(2,:) = box(2,:) + 400*eps('double');
        end
    else
        if size(box, 2) < dim
            disp 'Inconsistent box dimension with the given dim'
            parent = []; children = []; cluster = []; permutation = []; n = 0;
            return ;
        else
            %   minor modification to avoid rounding error
            box(1,:) = box(1,:) - 200*eps('double');
            box(2,:) = box(2,:) + 400*eps('double');
        end     
    end
    

    
    %%   Leaf node cluster with small number of points
    if (N <= minNum)
        parent = [0];
        children = NaN(1, nchildren);
        cluster = [1,N];
        enbox{1} = box;
        permutation = (1:N)';
        n = 1; 
        return ;
    end
    
    %%   Leaf node cluster with small box dimension
    idx_dim = find(box(2, :) > minBox);
    if (isempty(idx_dim)) % max(box(2,:)) < minBox
        parent = [0];
        children = NaN(1, nchildren);
        cluster = [1,N];
        enbox{1} = box;
        permutation = (1:N)';
        n = 1; 
        return ;
    end
    
    %%   Bisect the box along each long-enough dimension. 
    sub_clu = {};
    sub_box = {};
    rela_coord = bsxfun(@minus, coord, box(1, :));
    rela_index = floor(2 * bsxfun(@rdivide, rela_coord(:, idx_dim), box(2, idx_dim)));
    rela_index(rela_index == 2) = 1;    %  points on the box surface

    %   index of the sub-domain [i,j,k] is calculated  as 1+i*2^0+j*2^1+k*2^2
    idx = round(1 + sum(bsxfun(@times, rela_index, 2.^(0:length(idx_dim)-1)), 2));
    count = 0;
    for i = 1 : max(idx)
        tmp = find(idx == i);
        if isempty(tmp) 
            continue;
        else
            count = count + 1;
            sub_clu{count} = tmp; 
            tmpbox = box;
            tmpbox(:, idx_dim) = [box(1,idx_dim)+0.5*rela_index(tmp(1),:).*box(2,idx_dim)-20*eps('double'); 
                            0.5*box(2,idx_dim)+40*eps('double')];   
            sub_box{count} = tmpbox;         
        end                            
    end   
    
    %%   Recursively partition each subclusters obtained. 
    subb_par = cell(count,1);
    subb_chi = cell(count,1);
    subb_clu = cell(count,1);
    subb_box = cell(count,1);
    subb_per = cell(count,1);
    subb_n   = zeros(count,1);
    for i = 1 : count
        [subb_par{i}, subb_chi{i}, subb_clu{i}, subb_box{i}, subb_per{i}, subb_n(i)] = ...
            partition_bisec(coord(sub_clu{i},:), minSize, dim, sub_box{i});
    end
    
    
    %%   Combine the obtained children subtrees
    n = sum(subb_n) + 1;
    permutation = zeros(N,1);
    cluster = zeros(0, 2);
    idx_accu = 0;

    for i = 1 : count        
        permutation(sub_clu{i}) = subb_per{i} + idx_accu;        
        cluster = vertcat(cluster, subb_clu{i} + idx_accu);
        idx_accu = idx_accu + length(sub_clu{i});
    end
    cluster = [cluster; 1,N];    
    
    %%   Post-Order
    node_accu = 0;
    parent = zeros(0,1);
    children = zeros(0, nchildren);
    enbox = cell(0,1);
    for i = 1 : count
        parent = vertcat(parent, [subb_par{i}(1:end-1) + node_accu; n]);
        children = vertcat(children, subb_chi{i} + node_accu);
        node_accu = node_accu + subb_n(i);
        enbox = vertcat(enbox, subb_box{i});
    end
    parent = [parent; 0];
    children = [children; cumsum(subb_n'), NaN(1, nchildren-count)];
    enbox = vertcat(enbox, box); 
end

% % 
% % %   Non-leaf node with each dimension larger than the threhold 'minBox'
% %     if (length(idx_dim) == dim)        
% %         sub_clu = {};
% %         sub_box = {};
% %         rela_coord = bsxfun(@minus, coord, box(1,:));
% %         rela_index = floor(2 * bsxfun(@rdivide, rela_coord, box(2, :)));
% %         rela_index(rela_index == 2) = 1;    %   on the box surface
% % 
% %         %   index of the sub-domain [i,j,k] is calculated  as 1+i*2^0+j*2^1+k*2^2
% %         idx = 1 + sum(bsxfun(@times, rela_index, 2.^(0:dim-1)), 2);
% %         count = 0;
% %         for i = 1 : nchildren
% %             tmp = find(idx == i);
% %             if isempty(tmp) 
% %                 continue;
% %             else
% %                 count = count + 1;
% %                 sub_clu{count} = tmp; 
% %                 sub_box{count} = [box(1,:) + 0.5*rela_index(tmp(1),:).*box(2,:)-2*eps('double'); 
% %                                 .5*box(2,:)+4*eps('double')];            
% %             end                            
% %         end              
% %     else
% %         sub_clu = {};
% %         sub_box = {};
% %         rela_coord = bsxfun(@minus, coord, box(1,:));
% %         rela_index = floor(2 * bsxfun(@rdivide, rela_coord(:, idx_dim), box(2, idx_dim)));
% %         rela_index(rela_index == 2) = 1;    %   on the box surface
% % 
% %         %   index of the sub-domain [i,j,k] is calculated  as 1+i*2^0+j*2^1+k*2^2
% %         idx = 1 + sum(bsxfun(@times, rela_index, 2.^(0:length(idx_dim)-1)), 2);
% %         count = 0;
% %         for i = 1 : nchildren
% %             tmp = find(idx == i);
% %             if isempty(tmp) 
% %                 continue;
% %             else
% %                 count = count + 1;
% %                 sub_clu{count} = tmp; 
% %                 tmpbox = box;
% %                 tmpbox(:, idx_dim) = [box(1,idx_dim)+0.5*rela_index(tmp(1),:).*box(2,idx_dim)-2*eps('double'); 
% %                                 0.5*box(2,idx_dim)+4*eps('double')];   
% %                 sub_box{count} = tmpbox;         
% %             end                            
% %         end   
% %     end

