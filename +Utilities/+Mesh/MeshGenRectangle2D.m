function [EToV, VX, VY] = MeshGenRectangle2D(np, start_coor, end_coor)
% Input:
%   np - No of vertice on one edge
%   [a, b] - interval

m = np - 1; % elements on each edge

VX = linspace(start_coor, end_coor, np); 
VY = linspace(start_coor, end_coor, np)'; 

VX = repmat(VX, 1, np);
VY = repmat(VY, 1, np)'; 
VY = VY(:);

EToV = zeros(m^2, 4);
for irow = 1:m
    for icol = 1:m
        % Counterclockwise
        index = (irow-1)*m + icol;
        EToV(index, :) =  ...
            [np*(irow-1)+icol, np*(irow-1)+icol+1, ...
            np*irow+icol+1, np*irow + icol, ];
    end% for
end% for

end% func