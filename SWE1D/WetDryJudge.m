function isWet = WetDryJudge(mesh, h, physics)
% 

hDry = physics.getVal('minDepth');

hmean = CellMean(mesh, h);
% wetNode = (h > hDry); % bool value of wet nodes
% isWet = any(wetNode); % bool value of wet cells
isWet = hmean > hDry;
end% func