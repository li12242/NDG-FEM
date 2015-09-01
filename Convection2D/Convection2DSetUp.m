function var = Convection2DSetUp
% 2D convection problem
% dc/dt + udc/dx + vdu/dy = 0

N = 3;
% read triangle mesh
[EToV, VX, VY, EToR, BC] = Utilities.MeshReaderTriangle('Convection2D/mesh/rectangle');

tri = StdRegions.Triangle(N);
mesh = MultiRegions.RegionTriBC(tri, EToV, VX, VY, BC);
var = ConvectionInit(mesh);

Speed = [1,0]; % speed of domain, [u, v]
FinalTime = pi;

var = Convection2DSolver(mesh, var, FinalTime, Speed);

postprocess(mesh, var);
end% func

function postprocess(mesh, var)
plot3([mesh.x; mesh.x(1,:)], [mesh.y; mesh.y(1,:)],...
    [var; var(1,:)], '.')

% patch('Faces',EToV,...
%     'Vertices',[mesh.x(:), mesh.y(:)],...
%     'edgecol','k','facecol',[.8,.9,1])
end% func

function var = ConvectionInit(mesh)
% var = zeros(size(mesh.x));
var = sin(2*pi*mesh.x);
end% func