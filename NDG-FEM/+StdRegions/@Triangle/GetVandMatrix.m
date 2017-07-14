function Va = GetVandMatrix(obj, N, r, s)
% GETVANDMATRIX ����Vandermonde����
% Vandermonde����Ԫ������ V_{ij} = \psi_j(\xi_i)

Np = (N+1)*(N+2)/2;
Va = zeros(numel(r), Np);

for ind = 1:Np
    Va(:, ind) = OrthogonalFun(r,s,N,ind);
end