function flag = domain_box_hole(coord, box, hole)

dim = size(box, 1);
N = size(coord, 1);
flag_box = true(N, 1);
flag_hole = true(N, 1);
for i = 1 : dim
    flag_box = flag_box & bsxfun(@ge, coord(:, i), box(i, 1)) & bsxfun(@le, coord(:, i), box(i, 2));
    flag_hole = flag_hole & bsxfun(@ge, coord(:, i), hole(i, 1)) & bsxfun(@le, coord(:, i), hole(i, 2));
end
flag = flag_box & ~flag_hole;
end