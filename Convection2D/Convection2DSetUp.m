function var = Convection2DSetUp
% 2D convection problem
% dc/dt + udc/dx + vdu/dy = 0

N = 1;
mesh = quadSolver(N);
var = ConvectionInit(mesh);

Speed = [1,0]; % speed of domain, [u, v]
FinalTime = 1;

var = Convection2DSolver(mesh, var, FinalTime, Speed);

% postprocess(mesh, var, VX, VY, EToV);
plot3(mesh.x(mesh.vmapP),mesh.y(mesh.vmapP), var(mesh.vmapM))
end% func

function mesh = triSolver(N)

% read triangle mesh
[EToV, VX, VY, EToR, BC] = Utilities.Mesh.MeshReaderTriangle('Convection2D/mesh/triangle');

tri = StdRegions.Triangle(N);
mesh = MultiRegions.RegionTri(tri, EToV, VX, VY);

end% func

function mesh = quadSolver(N)

% read triangle mesh
[EToV, VX, VY, EToR, BC] = Utilities.Mesh.MeshReaderQuad('Convection2D/mesh/quad');

% uniform mesh
% [EToV, VX, VY] = uniformRetangle;

temp = EToV(:, 3); 
EToV(:, 3) = EToV(:, 4);
EToV(:, 4) = temp;

quad = StdRegions.Quad(N);
mesh = MultiRegions.RegionQuad(quad, EToV, VX, VY);

end% func

function [EToV, VX, VY] = uniformRetangle
np = 5; % nodes on each row
m = np - 1; % elements on each row

VX = 1:np; VX = repmat(VX, 1, np);
VY = [np:-1:1]'; VY = repmat(VY, 1, np)'; VY = VY(:);

EToV = zeros(m^2, 4);
for irow = 1:m
    for icol = 1:m
        index = (irow-1)*m + icol;
        EToV(index, :) = [np*(irow-1)+icol, np*(irow-1)+icol+1, ...
            np*(irow)+icol, np*(irow)+icol+1];
    end
end

end% func

function postprocess(mesh, var, VX, VY, EToV)

interp = TriScatteredInterp(mesh.x(:), mesh.y(:), var(:));
c = interp(VX, VY);

patch('Faces',EToV,...
    'Vertices',[VX, VY, c],...
    'edgecol','k','facecol',[.8,.9,1])
end% func

function var = ConvectionInit(mesh)
% var = ones(size(mesh.x));
var = mesh.x;
% xc = mean(mesh.x);
% left = xc < 0.5; right = xc > 0.5;
% var = sin(pi*mesh.x);%.*sin(2*pi*mesh.y);
% var = zeros(size(mesh.x));
% var(:,left) = 1; var(:,right) = 0;
end% func