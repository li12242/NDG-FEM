function refineflag = RefinedCellIdentify(mesh, h, bedElva)
% identify the wet/dey interface element
% Input: 
%   h - water depth
%   bedElva - bottom elevation
% Output:
%   transIndex - bool flag for wet/dry transition element, size [1, Ne]

% hPositive = 10^-3;

% identify transitation element
transIndex = TransiteCellIdentify(mesh, h);

% refinement condition - partically wet
deltaB = max( bedElva ) - min( bedElva );
hmean = CellMean(mesh,h);
refineflag = 2*hmean < deltaB;
% bedMean = CellMean(mesh, bedElva);
% bedMax = max(bedElva);
% refineflag = (hmean + bedMean + eps) < (bedMax);

refineflag = refineflag(:) & transIndex(:);
end% function 