function admnode = H2__admissible_nodelist(htree, near, far)

%   For each cluster i at kth level, find admissible node list containing 
%   node j if it satisfies either of the following conditions
%       1. j is at kth level and i*j is an admissible pair <= one pair of
%       ancestors of i and j are in the 'far' list. 
%       2. j is a leafnode at higher level (e.g., i's parent's level)
%       3. j is at lower level. This only happens when i is a
%       leafnode. In this case, (i,j) should be in the 'far' list. 
%   
%   For each cluster i, collect the index of the columns, colidx{i}, that
%   are admissible (being low-rank at the A(cluster(i), :) block rows). 

%   Tree info
parent   = htree.parent;
children = htree.children;
level    = htree.level;
nlevel   = length(level);


%   Collect admissible nodes information
admnode = cell(length(parent), 1);
for i = 2 : nlevel
   for j = 1 : length(level{i})
       node = level{i}(j);
       
       %    admissible block inherent from parent
       tmp = {};
       par_admclst = admnode{parent(node)};
       for k = 1 : length(par_admclst)
           if isnan(children(par_admclst(k),1)) %leafnode
               tmp{end+1} = par_admclst(k);
               continue;
           else         %non-leaf node
               nn = find(~isnan(children(par_admclst(k),:)), 1, 'last');
               tmp{end+1} = children(par_admclst(k),1:nn);
           end
       end       
       
       %    admissible pair at the same level
       %    ...if 'node' is a leafnode, this list also possibly contains
       %    nodes at lower-level, simply include those. 
       farnode = [far(far(:,2)==node, 1)', far(far(:,1)==node,2)']; 
       admnode{node} = horzcat(tmp{:}, farnode);             
   end
end