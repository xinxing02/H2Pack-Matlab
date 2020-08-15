function flag = domain_box(coord, box)

dim = size(box, 1);
N = size(coord, 1);
flag = true(N, 1);
for i = 1 : dim
    flag = flag & bsxfun(@ge, coord(:, i), box(i, 1)) & bsxfun(@le, coord(:, i), box(i, 2));
end
end