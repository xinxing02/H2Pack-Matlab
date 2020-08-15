function ulvfactor = HSS_ULV_LU(hssmat, htree, shift)
%
%   Construct the ULV decomposition of a general symmetric HSS matrix. 
%

%   Q:       Orthogonal matrix w.r.t. each node's basis
%   Lr, Lc:  LU factor w.r.t. each node's diagonal block
%   Idx:     Record the index of the submatrix which Q and L performs on for
%            each node in global sense. 

    %   diagonal shift
    if (nargin < 5)
        shift = 0;
    end
    
    %%  Basic information of the h2mat
    kdim = hssmat.kdim;
    if hssmat.JIT == true
        hssmat = H2__denseblocks(hssmat, htree);
    end
    D = hssmat.D;
    B = hssmat.B;
    U = hssmat.U;
    
    %%  Basic information of the htree
    children = htree.children;
    root     = htree.root;
    level    = htree.level;
    nlevel   = length(level);
    nnode    = htree.nnode;
    mcluster = [kdim * htree.cluster(:,1) - (kdim-1), kdim * htree.cluster(:,2)];

    %%  Initialization
    Q = cell(nnode, 1);         % Orthogonal matrix
    Lr = cell(nnode, 1);        % LU factor
    Lc = cell(nnode, 1);        % LU factor
    D_mid = cell(nnode, 1);     % Diagonal block 
    U_mid = cell(nnode, 1);     % Basis changed
    
    %%  Non-root node construction
    % In HSS, the top level with admissible blocks is the 2nd level. 
    for i = nlevel : -1 : 2         
        for j = 1 : length(level{i})
            node = level{i}(j); 
            child_node = children(node, ~isnan(children(node, :)));
            
            %   FIRST CASE : leaf node
            if isempty(child_node)  
                %   Size information of U_node and Diag_node
                nrow = size(D{node, node},1);
  
                %   QR decomposition of U{node}
                [Q{node}, tmpU] = qr(U{node});
                cutpoint = min(size(U{node}));
                
                %   U_mid{node} to be the triangular matrix only. 
                U_mid{node} = tmpU(1:cutpoint,:);
                
                %   Apply Q{node} on both sides of D{node, node}
                tmpD = Q{node}' * (D{node,node}+shift*eye(size(D{node,node}))) * Q{node};
                %   NOTE: tmpD is split into 2*2 blocks [tmpD11, tmpD12; tmpD21, tmpD22] with
                %   tmpD11 of dimension cutpoint*cutpoint.
                
                tmpD11 = tmpD(1:cutpoint, 1:cutpoint);
                tmpD12 = tmpD(1:cutpoint, (cutpoint+1):end);
                tmpD21 = tmpD12';
                tmpD22 = tmpD((cutpoint+1):end, (cutpoint+1):end);
                
                
                %   LU Decomposition of tmpD22               
                [tmpLr, tmpLc, tmpP] = lu(tmpD22);
                tmpLr = tmpP' * tmpLr; 
                %   UNDONE: tmpLr is no longer lower-triangular, need to be
                %   properly addressed in C. 
                
                %   singular case
                if (min(abs(diag(tmpLc))) < eps('double'))
                    disp('The target matrix is near-singular, consider adding a shift\n');
                    ulvfactor = [];
                    return 
                end
                
                %   Gaussian eliminate the off-diagonal blocks of tmpD, 
                %   i.e., eliminate tmpD12, and tmpD21 using tmpD22
                %   CASE 1: tmpD22 is not empty
                if size(U{node}, 2) < size(U{node}, 1) 
                    LD21 = tmpLr \ tmpD21;
                    DU12 = tmpD12 / tmpLc;
                %   CASE 2: tmpD22 is empty
                else
                    LD21 = zeros(0, nrow);
                    DU12 = zeros(nrow, 0);
                end        
                Lr{node} = [eye(cutpoint), DU12; zeros(nrow-cutpoint, cutpoint), tmpLr];
                Lc{node} = [eye(cutpoint), zeros(cutpoint, nrow-cutpoint); LD21, tmpLc];
                
                %   Schur complement at tmpD11.
                D_mid{node} = tmpD11 - DU12 * LD21;
                
                %   Row indices where Lr{node}, Lc{node} and Q{node} are applied to.
                Idx{node} = mcluster(node,1):mcluster(node,2);
            
            %   SECOND CASE : nonleaf node                    
            else 
                
                %   dimension of the each compressed children diagonal block.
                rowdim = arrayfun(@(x)size(D_mid{x},1), child_node);
                offset = [1, cumsum(rowdim)+1];
                nrow = sum(rowdim);
                
                %   build the compressed diagonal block
                tmpD = zeros(nrow);
                for k = 1 : length(child_node)
                    %   diagonal blocks
                    c1 = child_node(k);
                    idx1 = offset(k):offset(k+1)-1;
                    tmpD(idx1, idx1) = D_mid{c1};
                    %   off-diagonal blocks
                    for l = (k+1) : length(child_node)
                        c2 = child_node(l);
                        idx2 = offset(l):offset(l+1)-1;
                        tmpD(idx1, idx2) = U_mid{c1} * B{c1, c2} * U_mid{c2}';
                        tmpD(idx2, idx1) = tmpD(idx1, idx2)';
                    end
                end
                
                %   build the compressed basis matrix
                tmpU = blkdiag(U_mid{child_node}) * U{node};
                            
                %   QL decomposition of the compressed basis tmpU
                [Q{node}, tmpU] = qr(tmpU);
                cutpoint = min(size(tmpU));
 
                %   U_mid{node} to be the triangular matrix only. 
                U_mid{node} = tmpU(1:cutpoint,:);
                
                %   Apply Q{node} on both sides of tmpD
                tmpD = Q{node}' * tmpD * Q{node};
 
                tmpD11 = tmpD(1:cutpoint, 1:cutpoint);
                tmpD12 = tmpD(1:cutpoint, (cutpoint+1):end);
                tmpD21 = tmpD12';
                tmpD22 = tmpD((cutpoint+1):end, (cutpoint+1):end);
                
                
                %   LU Decomposition of tmpD22              
                [tmpLr, tmpLc, tmpP] = lu(tmpD22);
                tmpLr = tmpP' * tmpLr;
                %   UNDONE: tmpLr is no longer lower-triangular, need to be
                %   properly addressed in C. 
                
                %   singular case
                if (min(abs(diag(tmpLc))) < eps('double'))
                    disp('The target matrix is near-singular, consider adding a shift\n');
                    ulvfactor = [];
                    return 
                end
                
                %   Gaussian eliminate the off-diagonal blocks of tmpD, 
                %   i.e., eliminate tmpD12, and tmpD21 using tmpD22
                %   CASE 1: tmpD22 is not empty
                if size(tmpU, 2) < size(tmpU, 1) 
                    LD21 = tmpLr \ tmpD21;
                    DU12 = tmpD12 / tmpLc;
                %   CASE 2: tmpD22 is empty
                else
                    LD21 = zeros(0, nrow);
                    DU12 = zeros(nrow, 0);
                end        
                Lr{node} = [eye(cutpoint), DU12; zeros(nrow-cutpoint, cutpoint), tmpLr];
                Lc{node} = [eye(cutpoint), zeros(cutpoint, nrow-cutpoint); LD21, tmpLc];
                
                %   Schur complement at tmpD11.
                D_mid{node} = tmpD11 - DU12 * LD21;
 
                %   Row indices where Lr{node}, Lc{node} and Q{node} are applied to.
%                 tmpidx = arrayfun(@(x)Idx{x}(1:size(D_mid{x},1)), child_node, 'UniformOutput', false);
                tmpidx = cell(length(child_node), 1);
                for kk = 1 : length(child_node)
                    cnode = child_node(kk);
                    tmpidx{kk} = Idx{cnode}(1:size(D_mid{cnode},1));
                end
                Idx{node} = horzcat(tmpidx{:});                
            end
        end
    end

    %%   Last Step : root node
    child_node = children(root, ~isnan(children(root, :)));
    
    %   dimension of the each compressed children diagonal block.
    rowdim = arrayfun(@(x)size(D_mid{x},1), child_node);
    offset = [1, cumsum(rowdim)+1];
    nrow = sum(rowdim);

    %   build the compressed diagonal block
    D_mid{root} = zeros(nrow);
    for k = 1 : length(child_node)
        %   diagonal blocks
        c1 = child_node(k);
        idx1 = offset(k):offset(k+1)-1;
        D_mid{root}(idx1, idx1) = D_mid{c1};
        %   off-diagonal blocks
        for l = (k+1) : length(child_node)
            c2 = child_node(l);
            idx2 = offset(l):offset(l+1)-1;
            D_mid{root}(idx1, idx2) = U_mid{c1} * B{c1, c2} * U_mid{c2}';
            D_mid{root}(idx2, idx1) = D_mid{root}(idx1, idx2)';
        end
    end

    %   LU decomposition
    [Lr{root},Lc{root}] = lu(D_mid{root});    
    if (min(abs(diag(Lc{root}))) < eps('double'))
        disp('The target matrix is near-singular, consider adding a shift\n');
        ulvfactor = [];
        return 
    end
    
    %   Orthogonalization (only for simplicity)
    Q{root} = eye(size(D_mid{root}));
    
    %   Operating row indices
    tmpidx = arrayfun(@(x)Idx{x}(1:size(D_mid{x},1)), child_node, 'UniformOutput', false);
    Idx{root} = horzcat(tmpidx{:}); 
    
    
    %%  Wrapup
    ulvfactor.Q = Q;
    ulvfactor.Lr = Lr;
    ulvfactor.Lc = Lc;
    ulvfactor.Idx = Idx;
    ulvfactor.LU = true;
    ulvfactor.Chol = false;
    ulvfactor.shift = shift;
end