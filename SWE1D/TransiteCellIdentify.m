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
wetIndex = any(iswet); 
% when adjacent element possess different wet/dry status
% transIndex assignment is true
transIndex = xor(wetIndex(mesh.EToE(:, 1)), wetIndex(mesh.EToE(:, 2)));
% transitation element shoule be wet
transIndex = transIndex & wetIndex;

end% function 

