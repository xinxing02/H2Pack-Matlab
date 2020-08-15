function inadmnode = H2__inadmissible_nodelist(htree, near, far)

%   For each cluster i at kth level, find the inadmissible node list containing 
%   node j if it satisfies either of the following condition
%       1. if j is at kth level, then i * j must be inadmissible. 
%       2. j is at higher level. then j must be a leaf node and i * j must
%       be inadmissible => par(i) must be inadmissible to j too. 
%   
%   For each cluster i, collect the index of the columns, colidx{i}, that
%   are inadmissible (being low-rank at the A(cluster(i), :) block rows). 

%   Tree info
parent   = htree.parent;
level    = htree.level;
nodelvl  = htree.nodelvl;
minlvl   = 2;
nlevel   = length(level);
leafnode = htree.leafnode;
leaflvl  = nodelvl(leafnode);

admnode   = H2__admissible_nodelist(htree, near, far);
inadmnode = cell(length(parent), 1);

for i = nlevel : -1 : minlvl
    leaf2check = leafnode(leaflvl < i);
    candidate = [level{i}, leaf2check(:)'];
    for j = 1 : length(level{i})
        node = level{i}(j);
        inadmnode{node} = setdiff(candidate, [admnode{node},node]);
    end
end


