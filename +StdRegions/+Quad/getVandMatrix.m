function VandMatrix = getVandMatrix(nOrder, r, s)
% get Vandermonde matrix of quadrilateral element V
% $V_{ij} = \psi_j(\mathbf{\xi_i})$
% where $\psi_{(n+1)i+j+1}(\mathbf{\xi}) = P_i{r}P_j(s), \quad 0\le i, j \le n$
% is the orthogonal basis, and $P_i(x)$ is the normalized Legendre polynomials

VandMatrix = zeros(numel(r), (nOrder+1)^2 );
np = nOrder + 1; % number of points on edge
for i = 0:nOrder
    % P_{j-1}(r_i)$
    temp = Polylib.JacobiP(r(:), 0, 0, i);
    VandMatrix(:,i*np+1:(i+1)*np) = repmat(temp, 1, np);
    for j = 0:nOrder
        VandMatrix(:,i*np+1 + j) = VandMatrix(:,i*np+1 + j).*Polylib.JacobiP(s(:), 0, 0, j);
    end% for
end% for

end% func