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
hP = max( h(mesh.vmapP) );
hmean = CellMean(mesh,h);
refineflag = 2*hmean + hPositive < hP;
% bedMean = CellMean(mesh, bedElva);
% bedMax = max(bedElva);
% refineflag = (hmean + bedMean + eps) < (bedMax);

refineflag = refineflag(:) & transIndex(:);
end% function 

