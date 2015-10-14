function [EToV, VX, VY, BCType] = MeshReaderGmsh2D(fileName)
fig = fopen(fileName);

% read till '$Nodes'
temp = fgetl(fig);
while ~strncmp(temp, '$Nodes', 7)
    temp = fgetl(fig);
end
% read vertices
nVertice = fscanf(fig, '%d\n', 1);
node = fscanf(fig, '%d %f %f %f\n', [4, nVertice]);
node = node';
VX = node(:,2); VY = node(:,3);
% read till '$Elements'
while ~strncmp(temp, '$Elements', 8)
    temp = fgetl(fig);
end
% read element


%finish
fclose(fig);
end