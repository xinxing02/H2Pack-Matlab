function pos = gridpoint_on_boxsurface(corner, dir, L, nsample, dim) 
%
%   uniform sampling over the surface of a given cuboid 
%   Input: 
%       dim:     point dimension
%       corner:  dim * 1, the minimum corner of the box
%       dir:     dim * dim, the three directions that describes the cuboid
%       L:       dim * 1, the edge length along each direction
%       nsample: number of points to sample on the surface
%

if dim == 3
    Lx = L(1); Ly = L(2); Lz = L(3);
    corner = reshape(corner, 1, []);
    nsample_xy = ceil(0.5 * nsample * (Lx*Ly) / (Lx*Ly + Lx*Lz + Ly*Lz));
    Nx = ceil(sqrt(nsample_xy * Lx / Ly));
    Ny = ceil(Nx * Ly / Lx);
    Nz = ceil(Nx * Lz / Lx);
    nsample_xy = Nx * Ny;
    nsample_yz = Ny * Nz;
    nsample_xz = Nx * Nz;

    pos = zeros(2*nsample_xy+2*nsample_xz+2*nsample_yz, 3);
    idx = 0;
    intX = 1/(Nx + 1) * (1 : Nx)';
    intY = 1/(Ny + 1) * (1 : Ny)';
    intZ = 1/(Nz + 1) * (1 : Nz)';


    %   x-y plane
    pos(idx + (1:nsample_xy), :) = bsxfun(@plus, corner, ...
        Lx * kron(ones(Ny,1), intX) * dir(:,1)' + Ly * kron(intY, ones(Nx,1)) * dir(:,2)');
    idx = idx + nsample_xy;
    pos(idx + (1:nsample_xy), :) = bsxfun(@plus, corner + Lz*dir(:,3)', ...
        Lx * kron(ones(Ny,1), intX) * dir(:,1)' + Ly * kron(intY, ones(Nx,1)) * dir(:,2)');
    idx = idx + nsample_xy;

    %   x-z plane
    pos(idx + (1:nsample_xz), :) = bsxfun(@plus, corner, ...
        Lx * kron(ones(Nz,1), intX) * dir(:,1)' + Lz * kron(intZ, ones(Nx,1)) * dir(:,3)');
    idx = idx + nsample_xz;
    pos(idx + (1:nsample_xz), :) = bsxfun(@plus, corner + Ly*dir(:,2)', ...
        Lx * kron(ones(Nz,1), intX) * dir(:,1)' + Lz * kron(intZ, ones(Nx,1)) * dir(:,3)');
    idx = idx + nsample_xz;

    %   y-z plane
    pos(idx + (1:nsample_yz), :) = bsxfun(@plus, corner, ...
        Ly * kron(ones(Nz,1), intY) * dir(:,2)' + Lz * kron(intZ, ones(Ny,1)) * dir(:,3)');
    idx = idx + nsample_yz;
    pos(idx + (1:nsample_yz), :) = bsxfun(@plus, corner + Lx*dir(:,1)', ...
        Ly * kron(ones(Nz,1), intY) * dir(:,2)' + Lz * kron(intZ, ones(Ny,1)) * dir(:,3)');
end

if dim == 2
    Lx = L(1); Ly = L(2);
    corner = reshape(corner, 1, []);
    Nx = ceil(0.5 * nsample * Lx / (Lx + Ly));
    Ny = ceil(Nx * Ly / Lx);
    pos = zeros(2*(Nx + Ny), 2);
    idx = 0;
    intX = 1/(Nx + 1) * (1 : Nx)';
    intY = 1/(Ny + 1) * (1 : Ny)';
    
    %   x-direction
    pos(idx + (1:Nx), :) = bsxfun(@plus, corner, Lx * intX * dir(:,1)');
    idx = idx + Nx;
    pos(idx + (1:Nx), :) = bsxfun(@plus, corner + Ly * dir(:,2)', Lx * intX * dir(:,1)');
    idx = idx + Nx;
    
    %   y-direction
    pos(idx + (1:Ny), :) = bsxfun(@plus, corner, Ly * intY * dir(:,2)');
    idx = idx + Ny;
    pos(idx + (1:Ny), :) = bsxfun(@plus, corner + Lx * dir(:,1)', Ly * intY * dir(:,2)');
end


end
