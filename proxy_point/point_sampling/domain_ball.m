function flag = domain_ball(coord, center, radius)

dist = pdist2(center(:)', coord);
flag = (dist <= radius);
