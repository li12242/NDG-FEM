function isWet = SWE_WetEle1d(physics, mesh, h)
% 

hDry = physics.minDepth;

hmean = SWE_CellMean1d(mesh, h);
% wetNode = (h > hDry); % bool value of wet nodes
% isWet = any(wetNode); % bool value of wet cells
isWet = hmean > hDry;
end% func