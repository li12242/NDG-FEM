function VandMatrix=GetVandMatrix(~,nOrder,r)

% get Vandmonde Matrix of line element
% INPUT:
%   N   - highest order of base expansions P_N
%   r   - Lagrange node coordinate, between [-1,1]
% OUTPUT:
%   V   - Vandmonde Matrix, size: [size(r), N+1]
%         $V_{ij} = P_j(r_i)$
% USAGES: 
%   Vand = VandMatrix(N, r)

VandMatrix = zeros(numel(r),nOrder+1);
for j=0:nOrder
    % P_{j-1}(r_i)$
    VandMatrix(:,j+1) = Polylib.JacobiP(r(:), 0, 0, j);
end% for
end% func