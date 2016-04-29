function isWet = WetDryJudge(mesh, h, physics)
%

hDry = physics.getVal('minDepth');

wetNode = (h > hDry); % bool value of wet nodes
isWet = any(wetNode); % bool value of wet cells

end% func