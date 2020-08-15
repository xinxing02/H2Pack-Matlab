function flag = isadmissible_box(box1, box2, alpha)

%   checking the admissiblility of two boxes
%   use the L_\infty distance

%   Check each dimension separately
dim = size(box1, 2);
for i = 1 : dim
    if box1(2, i) < box2(2, i)  % in case of two boxes not at the same level. 
        dist = box1(1, i) - box2(1, i);
        if dist < 0 % box2 is on the right of box1 along ith dim
            if abs(dist) >= (1+alpha-1e-8) * box1(2, i)
                flag = true;
                return ;
            end
        else    % box2 is on the left of box1 along ith dim
            if dist >= (alpha-1e-8) * box1(2, i) + box2(2, i)
                flag = true;
                return ;
            end
        end
    else
        dist = box2(1, i) - box1(1, i);
        if dist < 0 %   box1 is on the right of box 2
            if abs(dist) >= (1 + alpha - 1e-8) * box2(2, i)
                flag = true;
                return ;
            end
        else %  box1 is on the leaf of box2
            if dist >= (alpha - 1e-8) * box2(2, i) + box1(2, i)
                flag = true;
                return ;
            end
        end
    end
end
flag = false;
return ;
end
