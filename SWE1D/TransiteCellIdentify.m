function transIndex = TransiteCellIdentify(mesh, isWet)
% identify the wet/dey interface element
% Input: 
%   mesh - mesh object
%   isWet - bool value of wet cells
% Output:
%   transIndex - bool flag for wet/dry transition element, size [1, Ne]

% % identify transitation element
% hPositive = 10^-3;
% % define wet cells
% iswet = (h > hPositive);
% isWet = any(iswet); 
% when adjacent element possess different wet/dry status
% transIndex assignment is true
transIndex = xor(isWet(mesh.EToE(:, 1)), isWet(mesh.EToE(:, 2)));
% transitation element shoule be wet
transIndex = transIndex & isWet;

end% function 

