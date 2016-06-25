function [EToV, VX, VY] = MeshGenRectangle2D(np, a, b)
% Input:
%   np - No of vertice on one edge
%   [a, b] - interval

m = np - 1; % elements on each edge

VX = linspace(a, b, np); VX = repmat(VX, 1, np);
VY = linspace(a, b, np)'; VY = repmat(VY, 1, np)'; 
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

% temp = EToV(:, 3); 
% EToV(:, 3) = EToV(:, 4);
% EToV(:, 4) = temp;
end% func