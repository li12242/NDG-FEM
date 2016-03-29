function refineflag = TransiteCellIdentify(mesh, h, bedElva)
% identify the wet/dey interface element
% Input: 
%   h - water depth
%   bedElva - bottom elevation
% Output:
%   refineflag - bool flag for wet/dry refinement element, size [1, Ne]

% identify transitation element
hPositive = 10^-3;
% define wet cells
iswet = (h > hPositive);
wetIndex = any(iswet); 
% when adjacent element possess different wet/dry status
% transIndex assignment is true
transIndex = xor(wetIndex(mesh.EToE(:, 1)), wetIndex(mesh.EToE(:, 2)));
% transitation element shoule be wet
transIndex = transIndex & wetIndex;

% refinement condition
hmean = CellMean(mesh,h);
bedMean = CellMean(mesh, bedElva);
bedMax = max(bedElva);
refineflag = (hmean + bedMean) < (bedMax + hPositive);

refineflag = refineflag(:) & transIndex(:);
end% function 

function hmean = CellMean(mesh, h)
% get mean depth in each cell
Np = mesh.Shape.nNode;
uh = mesh.Shape.VandMatrix\h; uh(2:Np,:)=0;
uavg = mesh.Shape.VandMatrix*uh; hmean = uavg(1,:);
end% func