function [time, position] = FrontTrack
% get the wet/dry front position

% output file
filename = 'SWE1DTsunamiRunup.nc';

% rebuild element variables
x1 = -500; x2 = 50000;
np = 2; nele = 1000;
[~, VX, ~, EToV] = Utilities.Mesh.MeshGen1D(x1, x2, nele);
line = StdRegions.Line(np-1);
mesh = MultiRegions.RegionLine(line, EToV, VX);


% get results
time = ncread(filename, 'time'); x = ncread(filename, 'x');
itime = 1;
h = ncread(filename, 'h', [1, itime],[inf, 1]);


position = zeros(size(time));
transIndex = TransiteCellIdentify(mesh, h);
position(itime) = getFrontPos(mesh, h, transIndex);

% % bottom topography
% bedElva = 5000 - 0.1*x;

for itime = 1:numel(time)
    h = ncread(filename, 'h', [1, itime],[inf, 1]);

    transIndex = TransiteCellIdentify(mesh, h);
    position(itime) = getFrontPos(mesh, h, transIndex);

%     if (numel(transIndex) > 1) | (isempty(transIndex))
%         error('transition element number error.')
%     end% if
    fprintf('Processing: %f ...\n', itime/numel(time))
end% for
end% func

function xw = getFrontPos(mesh, h, index)
nodeIndex = [index*2-1; index*2];
vx = mesh.x(mesh.vmapM(:, index)); dx = vx(2) - vx(1);
hmean = CellMean( mesh, h(nodeIndex) );
hP = h(mesh.vmapP(:, index ));

deltax = 2*(hmean*dx) ./(hP(2));
xw = vx(1) + deltax;
end

function transIndex = TransiteCellIdentify(mesh, h)
% identify the wet/dey interface element
% Input: 
%   h - water depth
%   bedElva - bottom elevation
% Output:
%   transIndex - bool flag for wet/dry transition element, size [1, Ne]

% identify transitation element
hPositive = 10^-3;
% define wet cells
iswet = (h > hPositive);
wetIndex = (iswet(1:2:end-1) + iswet(2:2:end)) > 0;
% when adjacent element possess different wet/dry status
% transIndex assignment is true
transIndex = xor(wetIndex(mesh.EToE(:, 1)), wetIndex(mesh.EToE(:, 2)));
% transitation element shoule be wet
transIndex = find(transIndex & wetIndex);
end% function 