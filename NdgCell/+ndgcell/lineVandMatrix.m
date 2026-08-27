function V = lineVandMatrix( N, r )
%> @brief 1D orthonormal Jacobi (Legendre) Vandermonde matrix.
%> @details \f$ V_{i,j+1} = \tilde P_j(r_i) \f$, \f$ j=0,\cdots,N \f$.
%> Shared by StdTri/StdPrismTri node warping.
% ======================================================================
%> This function is part of the NDG-FEM software.
% ======================================================================
V = zeros(numel(r), N+1);
for j = 0:N
    V(:, j+1) = JacobiP(r(:), 0, 0, j);
end% for
end% func
