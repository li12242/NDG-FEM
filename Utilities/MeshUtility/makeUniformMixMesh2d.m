function [ mesh ] = makeUniformMixMesh2d( N, xlim, ylim, Mx, My, cellType, bcType )

checkInput(xlim, ylim, Mx, My, bcType);
tri = StdTri(N);
quad = StdQuad(N);

flag = 0;
% Parameters
Nx = Mx + 1; % number of nodes along x coordinate
Ny = My + 1;
K  = Mx * My * 2;
Nv = Nx * Ny;
EToR = NdgRegionType.Normal * ones(K, 1, 'int8');

xmin = min(xlim); xmax = max(xlim);
ymin = min(ylim); ymax = max(ylim);
VX   = linspace(xmin, xmax, Nx) ;
VY   = linspace(ymin, ymax, Ny)';
VX   = repmat(VX, 1, Ny) ;
VY   = repmat(VY, 1, Nx)';
vx   = VX(:);
vy   = VY(:);

triEToV = zeros(4, Mx*My);
ind = 1;
for i = 1:My
    for j = 1:Mx
        % vertex index
        v1 = Nx*(i-1) + j;
        v2 = Nx*(i-1) + j + 1;
        v3 = Nx*i + j;
        v4 = Nx*i + j + 1;
        % Counterclockwise
        EToV(:, ind)=[v1, v2, v4, v3]';
        ind = ind +1;
    end% for
end% for


end

