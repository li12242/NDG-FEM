function testing_RegionTri


%% testing RegionTri
nOrder = 3;
tri = StdRegions.Triangle(nOrder);

% testing pyhsics mapping
% warning: the vertice of element in EToV must be counterclockwise
EToV = [3,1,2];
VX = [0, 2, 1]'; VY = [0, 0, sqrt(3)]';
mesh = MultiRegions.RegionTri(tri, EToV,VX,VY);
if any( abs(mesh.sJ-1.0)>10^-10 )
    error('edge mapping error')
else
    fprintf('edge mapping pass\n')
end
eleArea = sqrt(3);
stdEleArea = 2;
if any( abs(mesh.J - eleArea/stdEleArea)>10^-10 )
    error('element area mapping error')
else
    fprintf('element area mapping pass\n')
end

% testing facelist mapping
% 1   3   5
% |--/|--/
% | / | /
% |/--|/
% 2   4
EToV = [3,1,2; 3,2,4; 5,3,4];
VX = [0, 0, 1, 1, 2]'; VY = [1, 0, 1, 0,1]';
mesh = MultiRegions.RegionTri(tri, EToV,VX,VY);

% check vmapP vmapM
vn = [1,2,3,1];
figure; hold on;
for ie = 1:mesh.nElement
    plot(VX(EToV(ie, vn)), VY(EToV(ie, vn)), 'k')
end
plot(mesh.x, mesh.y,'ro')
plot(mesh.x(mesh.vmapM), mesh.y(mesh.vmapM), 'bo')
axis equal


%% testing RegionTriBC
nOrder = 3;
[EToV, VX, VY, ~, BC] = Utilities.MeshReaderTriangle('SWE2D/grid/untitled');

tri = StdRegions.Triangle(nOrder);
mesh = MultiRegions.RegionTriBC(tri, EToV, VX, VY, BC);


end