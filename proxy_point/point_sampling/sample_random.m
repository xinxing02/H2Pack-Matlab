function Y = sample_random(domain, N, box)

%   domain: function handle, 
%           input: points, 
%           output: positive number (in the domain)
%                   negative number (outside of the domain)
%   N: number of points
%   box: a bounding box domain (dim*2 matrix) that encloses the domain

dim = size(box, 1);
edge = box(:, 2) - box(:, 1);
Y = zeros(N, dim);
index = 0;
while 1
    Y_initial = bsxfun(@times, rand(2*N, dim), edge');
    Y_initial = bsxfun(@plus, Y_initial, box(:,1)');
    result = Y_initial(domain(Y_initial) > 0, :);
    if size(result, 1) + index >= N
        Y( (index+1) : N, :) = result(1 : (N-index), :);
        return ;
    else
        Y( (index+1) : (index + size(result, 1)), :) = result;
        index = index + size(result, 1);
    end    
end

end
