function testing_RegionLine

nOrder = 8;
line = StdRegions.Line(nOrder);

EToV = [1,2; 2,3; 3,5; 5,4];
VX = [1,2,3,5,4];

%% RegionLine
mesh = MultiRegions.RegionLine(line, EToV, VX);

% check points
y = (mesh.x).*0;
figure
plot(mesh.x,y,'or');

%% RegionLineBC
BC = [2, 1; 3, 4];
meshBC = MultiRegions.RegionLineBC(line, EToV, VX, BC);

% check points
y = (mesh.x).*0;
figure; hold on
plot(mesh.x,y,'or');
plot(mesh.x(mesh.vmapM(meshBC.mapI)), y(1), 'bo');
plot(mesh.x(mesh.vmapM(meshBC.mapO)), y(1), 'co');

end