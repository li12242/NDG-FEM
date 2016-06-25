function [EToV, VX, VY, EToR, BC] = MeshReaderQuad(fileHead)
% read quad element file
% Input:    fileHead - head name of file, .e.g 'SWE2D' for 
%                   'SWE2D.edge', 'SWE2D.ele' & 'SWE2D.node'
% Output:   EToV - element to vertice, size [nEle x 3]
%           EToR - element to regions, size [nEle x 1]
%           VX   - coordinate x
%           VY   - coordinate y
%           BC   - boundary condition, [BCType, node1, node2]

% read *.node file
nodefile = [fileHead, '.node'];
fig = fopen(nodefile);
temp = fscanf(fig, '%d %d %d %d\n', [4,1]);
nNode = temp(1); nDim = temp(2);
% get coordinate
Coor = fscanf(fig, '%d %f %f %f\n', [4,nNode]); Coor = Coor';
VX = Coor(:,2); VY = Coor(:,3);
fclose(fig);

% read *.ele file
elefile = [fileHead, '.ele'];
fig = fopen(elefile);
temp = fscanf(fig, '%d %d %d\n', [3,1]);
nElement = temp(1); nVerPerEle = temp(2);
Ele = fscanf(fig, '%d %d %d %d %d %d\n', [6,nElement]);
Ele = Ele';
EToV = Ele(:,2:5);
EToR = Ele(:,6);
fclose(fig);

% read *.edge file
edgefile = [fileHead, '.edge'];
fig = fopen(edgefile);
temp = fscanf(fig, '%d %d\n', [2,1]);
nBC = temp(1);
BC = fscanf(fig, '%d %d %d %d', [4, nBC]);
BC = BC'; BC = BC(:,[4,2,3]);
fclose(fig);
end