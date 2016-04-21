function [new_mesh, h1, q1, new_bedElva, localEleIndex] = ...
    Hrefine1D(mesh, refineflag, h, q, bedElva, physics)
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

%% get new vertex index and its location, EToV, h1, q1
Nv = length(VX(:)); 
Nrefine = sum(refineflag); % No. of refined elements

v1 = EToV(refineflag, 1); v2 = EToV(refineflag, 2); 
v3 = Nv + (1:Nrefine)'; % index of new vertex

% incert element to EToV
EToV(refineflag,:) = [v1, v3]; % left part
EToV((Ne+1):(Ne + Nrefine), :) = [v3, v2]; % right part

[new_mesh, newVX, newBedElva, h1, q1] = getRefinedValue(mesh, h, q, bedElva, refineflag);
VX(v3) = newVX;

refinedIndex = find(refineflag); % refined local element index
localEleIndex = [refinedIndex, [(Ne+1):(Ne + Nrefine)]' ];

%% get new EToE, EToF and other geometric variables
[new_mesh.EToE,new_mesh.EToF]= tiConnect1D(new_mesh, EToV);

[new_mesh.x, new_mesh.rx, new_mesh.J] = new_mesh.Shape.getEleGeometric(VX(EToV'));
[new_mesh.nx, new_mesh.sJ] = new_mesh.Shape.getFaceGeometric(new_mesh.x);

% [new_mesh.vmapM, new_mesh.vmapP] = new_mesh.BuildMap(EToV, new_mesh.EToE, new_mesh.EToF);
new_mesh = BuildMap(new_mesh, mesh, localEleIndex);

new_mesh.fScale = new_mesh.sJ./new_mesh.J(new_mesh.vmapM);
%% get new bottom elevation
VB = physics.getVal('VB');
VB(v3) = newBedElva;
vb = VB(EToV');
new_bedElva = 0.5*((1-mesh.Shape.r)*vb(1,:) + (mesh.Shape.r+1)*vb(2,:));

end% function

function new_mesh = BuildMap(new_mesh, mesh, localEleIndex)
ne = new_mesh.nElement; % # of elements
nfp = 1; % # of points on faces
nf = 2; % # of faces
np = new_mesh.Shape.nNode;
nrefine = size(localEleIndex, 1);

new_mesh.vmapM = zeros(nf*nfp, ne); new_mesh.vmapP = zeros(nf*nfp, ne);
new_mesh.vmapM(:, 1:mesh.nElement) = mesh.vmapM;

for ie = mesh.nElement+1:ne
    new_mesh.vmapM(:, ie) = (ie - 1)*np + new_mesh.Shape.getFaceListToNodeList;
end

new_mesh.vmapP(:, 1:mesh.nElement) = mesh.vmapP;
for ie = 1:nrefine
    k1 = localEleIndex(ie, 1); k2 = localEleIndex(ie, 2);
%     new_mesh.vmapP(1, k1) is the same with mesh
    new_mesh.vmapP(2, k1) = new_mesh.vmapM(1, k2);
    new_mesh.vmapP(1, k2) = new_mesh.vmapM(2, k1);
    new_mesh.vmapP(2, k2) = mesh.vmapP(2, k1);
    
    e2 = mesh.EToE(k1, 2); f2 = mesh.EToF(k1, 2); % f2 should be 1 (left side)
    new_mesh.vmapP(f2, e2) = new_mesh.vmapM(2, k2);
end% for

end% func

function [new_mesh, newVX, newBedElva, h1, q1] = ...
    getRefinedValue(mesh, h, q, bedElva, refineflag)
% calculate new vertex location
new_mesh = mesh;

Nrefine = sum(refineflag); % No. of refined elements
newVX = zeros(Nrefine, 1);
newBedElva = zeros(Nrefine, 1);
% allocate refined variables
h1 = zeros(size(h,1), Nrefine + mesh.nElement);
q1 = zeros(size(h1,1));
h1(:, 1:mesh.nElement) = h; q1(:, 1:mesh.nElement) = q;
% new_mesh.x = zeros(size(h1));
% new_mesh.x(:, 1:mesh.nElement) = mesh.x;
new_mesh.nElement = mesh.nElement + Nrefine;
% refined element index
eleIndex = find(refineflag);

for i = 1:Nrefine
    hP = h(mesh.vmapP([1, end], eleIndex(i) ));
    vx = mesh.x(mesh.vmapM([1, end], eleIndex(i) ));
    b = bedElva( mesh.vmapM([1, end], eleIndex(i)) );
    db = max(b) - min(b);
    hmean = CellMean(mesh, h(:, eleIndex(i) ));
    qmean = CellMean(mesh, q(:, eleIndex(i) ));
    dx = vx(2) - vx(1);
    
    if (hP(1) > hP(2)) && (hP(1) > 0) % left is wet cell
        deltax = 2*(hmean*dx) ./(hP(1));
        newVX(i) = vx(1) + deltax;
        newBedElva(i) = b(1) + (b(2) - b(1))/dx*deltax;
        % interpolate refined element value
        h1(:, eleIndex(i)) = 0.5*((1-mesh.Shape.r)*hP(1) ...
            + (mesh.Shape.r+1)*0);
        h1(:, mesh.nElement + i) = 0;
        q1(:, eleIndex(i)) = 0.5*((1-mesh.Shape.r)*2*qmean*dx/deltax);
        q1(:, mesh.nElement + i) = 0;
        
    elseif (hP(1) < hP(2)) && (hP(2) > 0) % right is wet
        deltax = 2*(hmean*dx) ./(hP(2)+0);
        newVX(i) = vx(2) - deltax;
        newBedElva(i) = b(2) - (b(2) - b(1))/dx*deltax;
        % interpolate refined element value
        h1(:, eleIndex(i)) = 0;
        h1(:, mesh.nElement + i) = 0.5*((1-mesh.Shape.r)*0 ...
            + (mesh.Shape.r+1)*hP(2));
        q1(:, eleIndex(i)) = 0;
        q1(:, mesh.nElement + i) = 0.5*((mesh.Shape.r+1)*2*qmean*dx/deltax);
    else % assum steady state
        if b(1) < b(2)
            vs = 0.5*(dx*db); vm = hmean*dx;
            deltax = dx*sqrt(vm/vs);
            newVX(i) = vx(1) + deltax;
            newBedElva(i) = b(1) + (b(2) - b(1))/dx*deltax;
            % interpolate refined element value
            h1(:, eleIndex(i)) = 0.5*((1-mesh.Shape.r)*2*vm./deltax ...
                + (mesh.Shape.r+1)*0);
            h1(:, mesh.nElement + i) = 0;
            q1(:, eleIndex(i)) = 0.5*((1-mesh.Shape.r)*2*qmean*dx/deltax);
            q1(:, mesh.nElement + i) = 0;
        elseif b(1) > b(2)
            vs = 0.5*(dx*db); vm = hmean*dx;
            deltax = dx*sqrt(vm/vs);
            newVX(i) = vx(2) - deltax;
            h1(:, eleIndex(i)) = 0;
            h1(:, mesh.nElement + i) = 0.5*((1-mesh.Shape.r)*0 ...
                + (mesh.Shape.r+1)*2*vm./deltax);
            q1(:, eleIndex(i)) = 0;
            q1(:, mesh.nElement + i) = 0.5*((mesh.Shape.r+1)*2*qmean*dx/deltax);
        else
            error('bottom is flat');
        end% if
    end% if
end% for
end% function

function hmean = CellMean(mesh, h)
% get mean depth in each cell
Np = mesh.Shape.nNode;
uh = mesh.Shape.VandMatrix\h; uh(2:Np,:)=0;
uavg = mesh.Shape.VandMatrix*uh; hmean = uavg(1,:);
end% func