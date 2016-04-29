function refineflag = RefinedCellIdentify(mesh, h, physics, isWet)
% identify the wet/dey interface element
% Input: 
%   h - water depth
%   bedElva - bottom elevation
% Output:
%   transIndex - bool flag for wet/dry transition element, size [1, Ne]

hPositive = physics.getVal('minDepth');

% identify transitation element
transIndex = TransiteCellIdentify(mesh, isWet);

% refinement condition - partically wet
hmean = CellMean(mesh,h);
hPmax = max(h(mesh.vmapP));
refineflag = ( hPmax > 0 ) & ( hPmax > 2*hmean ) & ( hmean >  hPositive);

refineflag = refineflag(:) & transIndex(:);
end% function 