function Vand = VandMatrix(N, r)
% Calculate 1D Vandmonde Matrix
% INPUT:
%   N   - highest order of base expansions P_N
%   r   - Lagrange node coordinate, between [-1,1]
% OUTPUT:
%   V   - Vandmonde Matrix, size: [size(r), N+1]
%         $V_{ij} = P_j(r_i)$
% USAGES: 
%   Vand = VandMatrix(N, r)
Vand = zeros(numel(r), N+1);
for j=0:N
    % P_{j-1}(r_i)$
    Vand(:,j+1) = Polylib.JacobiP(r(:), 0, 0, j);
end% for
end% func