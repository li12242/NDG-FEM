function [new_mesh, h1, q1, localEleIndex] = Hrefine1D(mesh, refineflag, h, q, physics)
% refine element
% Input:
%   refineflag - bool variable for defining the refined element, size [Ne, 1]
%   mesh - mesh object
%   h, q - variables
%   EToV - element to vertex
% Output:
%   new_mesh - new mesh object
%   h - refined variable of water depth
%   q - refined variable of discharge
%   localEleIndex - refined element index, size [No. of refined Element,2]

EToV = physics.getVal('EToV');
VX = physics.getVal('VX');
Ne = mesh.nElement; 
hDelta = 1e-3;

%% get new vertex index and its location, EToV, h1, q1
Nv = length(VX(:)); 
Nrefine = sum(refineflag); % No. of refined elements

v1 = EToV(refineflag, 1); v2 = EToV(refineflag, 1); 
v3 = Nv + (1:Nrefine)'; % index of new vertex

% incert element to EToV
EToV(refineflag,:) = [v1, v3]; % left part
EToV((Ne+1):(Ne + Nrefine), :) = [v3, v2]; % right part 

[new_mesh, newVX, h1, q1] = getRefinedValue(mesh, h, q, refineflag, hDelta);
VX(v3) = newVX;

%% get new EToE and EToF



end% function

function [new_mesh, newVX, h1, q1] = getRefinedValue(mesh, h, q, refineflag, hDelta)
% calculate new vertex location
new_mesh = mesh;

Nrefine = sum(refineflag); % No. of refined elements
newVX = zeros(Nrefine, 1);
% allocate refined variables
h1 = zeros(size(h,1), Nrefine + mesh.nElement);
q1 = zeros(size(h1,1));
h1(:, 1:mesh.nElement) = h; q1(:, 1:mesh.nElement) = q;
new_mesh.x = zeros(size(h1));
new_mesh.nElement = mesh.nElement + Nrefine;
% refined element index
eleIndex = find(refineflag);

for i = 1:Nrefine
    hP = h(mesh.vmapP(:, eleIndex(i) ));
    vx = mesh.x(mesh.vmapM(:, eleIndex(i) ));
    
    hmean = CellMean(mesh, h(:, eleIndex(i) ));
    qmean = CellMean(mesh, q(:, eleIndex(i) ));
    dx = vx(2) - vx(1);
    
    % element nodes location
    new_mesh.x(:, eleIndex(i)) = 0.5*((1-mesh.Shape.r)*vx(1) ...
        + (mesh.Shape.r+1)*newVX(i));
    new_mesh.x(:, mesh.nElement + i) = 0.5*((1-mesh.Shape.r)*newVX(i) ...
        + (mesh.Shape.r+1)*vx(2));
    
    if hP(1) > hP(2) % left is wet cell
        deltax = 2*(hmean*dx) ./hP(1);
        newVX(i) = vx(1) + deltax;
        % interpolate refined element value
        h1(:, eleIndex(i)) = 0.5*((1-mesh.Shape.r)*hP(1) ...
            + (mesh.Shape.r+1)*hDelta);
        h1(:, mesh.nElement + i) = hDelta;
        q1(:, eleIndex(i)) = 0.5*((1-mesh.Shape.r)*2*qmean);
        q1(:, mesh.nElement + i) = 0;
        
    else % right is wet
        deltax = 2*(hmean*dx) ./hP(2);
        newVX(i) = vx(2) - deltax;
        
        h1(:, eleIndex) = hDelta;
        h1(:, mesh.nElement + i) = 0.5*((1-mesh.Shape.r)*hDelta ...
            + (mesh.Shape.r+1)*hP(2));
        q1(:, eleIndex) = 0;
        q1(:, mesh.nElement + i) = 0.5*((mesh.Shape.r+1)*2*qmean);
        
    end% if
    
end% for
end% function

function hmean = CellMean(mesh, h)
% get mean depth in each cell
Np = mesh.Shape.nNode;
uh = mesh.Shape.VandMatrix\h; uh(2:Np,:)=0;
uavg = mesh.Shape.VandMatrix*uh; hmean = uavg(1,:);
end% func