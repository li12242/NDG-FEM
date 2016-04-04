function [h, q] = Hcombine1D(localEleIndex, h1, q1, refinedMesh, refineflag)
% combine local refined mesh
% Input: 
%   localEleIndex = refined element index, size [No. of refined Element,2]
% 
% 

Nrefine = size(localEleIndex, 1);

h = h1(:, 1:(end - Nrefine)); 
q = q1(:, 1:(end - Nrefine));

for i = 1:Nrefine
    e1 = localEleIndex(i, 1); e2 = localEleIndex(i, 2);
    eindex = min(e1, e2);
    
    % get left and right boundary value
    id = [refinedMesh.x(refinedMesh.vmapM(:, e1)), h1(refinedMesh.vmapM(:, e1)), ...
        q1(refinedMesh.vmapM(:, e1));...
        refinedMesh.x(refinedMesh.vmapM(:, e2)), h1(refinedMesh.vmapM(:, e2)), ...
        q1(refinedMesh.vmapM(:, e2))];
    id = sortrows(id);
    dx = id(4, 1) - id(1,1);
    hl = id(1, 2); %ql = id(1, 3);
    hr = id(4, 2); %qr = id(4, 3);
    % reinterp the conbine element value
    hmean1 = CellMean(refinedMesh, h1(:, e1));
    hmean2 = CellMean(refinedMesh, h1(:, e2));
    qmean1 = CellMean(refinedMesh, q1(:, e1));
    qmean2 = CellMean(refinedMesh, q1(:, e2));
    
    dx1 = refinedMesh.x(end, e1) - refinedMesh.x(1, e1);
    dx2 = refinedMesh.x(end, e2) - refinedMesh.x(1, e2);
    hmean = (hmean1*dx1 + hmean2*dx2)./dx;
    qmean = (qmean1*dx1 + qmean2*dx2)./dx;
   
    if hl > hr % left is wet
        h(:, eindex) = 0.5*( (1-refinedMesh.Shape.r)*2*hmean + (refinedMesh.Shape.r+1)*0 );
        q(:, eindex) = 0.5*( (1-refinedMesh.Shape.r)*2*qmean + (refinedMesh.Shape.r+1)*0 );
    else % right is wet
        h(:, eindex) = 0.5*( (1-refinedMesh.Shape.r)*0 + (refinedMesh.Shape.r+1)*2*hmean );
        q(:, eindex) = 0.5*( (1-refinedMesh.Shape.r)*0 + (refinedMesh.Shape.r+1)*2*qmean );
    end
    
end% for

h(:, ~refineflag) = h1(:, ~refineflag);
q(:, ~refineflag) = q1(:, ~refineflag);
end% func

function hmean = CellMean(mesh, h)
% get mean depth in each cell
Np = mesh.Shape.nNode;
uh = mesh.Shape.VandMatrix\h; uh(2:Np,:)=0;
uavg = mesh.Shape.VandMatrix*uh; hmean = uavg(1,:);
end% func