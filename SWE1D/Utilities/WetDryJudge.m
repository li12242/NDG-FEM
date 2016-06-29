%% WetDryJudge
% Determine wet/dry status
function wetEleFlag = WetDryJudge(phys, mesh, h)
% Input:
%   phys  - structure variable
%   mesh  - mesh object
%   h     - water depth
% Output:
%   wetEleFlag  - element flag for wet elements
% 
minDepth   = phys.minDepth;
hmean      = CellMean(mesh, h);
wetEleFlag = hmean > minDepth;
end% func