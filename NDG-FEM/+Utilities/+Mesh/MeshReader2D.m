function [EToV, VX, VY, EToR, BC] = MeshReader2D(fileHead)
% Read 2d mesh files.
% 
% Input:    fileHead - head name of file, .e.g 'SWE2D' for 
%                   'SWE2D.edge', 'SWE2D.ele' & 'SWE2D.node'
% Output:   EToV - element to vertice, size [nEle x 3]
%           EToR - element to regions, size [nEle x 1]
%           VX   - coordinate x
%           VY   - coordinate y
%           BC   - boundary condition, [BCType, node1, node2]
% Usages:
% 
% 

%% read *.node file
nodefile = [fileHead, '.node'];
fig   = fopen(nodefile);
temp  = fscanf(fig, '%d %d %d %d\n', [4,1]);
nNode = temp(1); % No. of vertex
nDim  = temp(2); % No. of dimension
if nDim ~= 2
    error('The dim of mesh file (%d) is incorrect', nDim)
end% if 
% get coordinate
Coor = fscanf(fig, '%d %f %f\n', [3,nNode]); 
VX   = Coor(2,:)'; 
VY   = Coor(3,:)';
fclose(fig);

%% read *.ele file
elefile = [fileHead, '.ele'];
fig     = fopen(elefile);
temp    = fscanf(fig, '%d %d %d\n', [3,1]);
Ne      = temp(1); % No. of element
Nvertex = temp(2); % No. of vertex in each element
Ele  = fscanf(fig, '%d %d %d %d %d\n', [Nvertex+2, Ne]);
Ele  = Ele';
EToV = Ele(:,(2:Nvertex+1));
EToR = Ele(:,Nvertex+2); % region index of each element
fclose(fig);

EToV = Utilities.Mesh.ResortVertex_Mex(EToV, VX, VY);

%% read *.edge file
edgefile = [fileHead, '.edge'];
fig      = fopen(edgefile);
temp     = fscanf(fig, '%d %d\n', [2,1]);
nBC      = temp(1);
BC = fscanf(fig, '%d %d %d %d', [4, nBC]);
BC = BC'; 
BC = BC(:,[4,2,3]);
fclose(fig);
end