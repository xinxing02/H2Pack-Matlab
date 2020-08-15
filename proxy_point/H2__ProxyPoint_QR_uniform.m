function Yp = H2__ProxyPoint_QR_uniform(kernel, htree, alpha, tol)
       
    %   Basic info
    level  = htree.level;
    enbox  = htree.enbox;
    nlevel = length(level);
    Yp     = cell(nlevel, 1);
    L = enbox{htree.root}(2, :);
    dim = length(L);
    
    %   selection of proxy points at each level
    maxL = enbox{level{1}(1)}(2, :)';
    for i = 3 : nlevel 
        %   properly initalize the two domains
        semi_L1 = enbox{level{i}(1)}(2, :)' / 2;
        box1 = [-semi_L1, semi_L1]; 
        hole2 = (1+2*alpha) * box1;
        semi_L3 = min([(1+8*alpha)*semi_L1, maxL - semi_L1], [], 2);
        box2 = [-semi_L3, semi_L3];
        domain1 = @(coord)(domain_box(coord, box1));
        domain2 = @(coord)(domain_box_hole(coord, box2, hole2));

        %   select the proxy points
        tmp_Yp = pp_numerical_selection_uniform(kernel, domain1, box1, domain2, box2, tol);
        if isempty(tmp_Yp)
            Yp{i} = zeros(0, dim);
        else
            Yp{i} = tmp_Yp;
            %   Can consider density if the accuracy is not good
            %   Yp{i} = sample_densify(tmp_Yp, 1, domain2);  
        end
    end
