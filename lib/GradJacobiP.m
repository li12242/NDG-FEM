function dP = GradJacobiP( r, alpha, beta, n )
%> @brief Derivative of the orthonormal Jacobi polynomial, \f$ d\tilde{P}_n^{(\alpha,\beta)}/dr \f$.
%>
%> Replaces the Polylib mex gateway (thirdParty/Polylib/GradJacobiP.c).
%> Uses the derivative identity of the normalized polynomial
%> (Hesthaven & Warburton, Eq. A.9):
%> \f[ \frac{d}{dr}\tilde{P}_n^{(\alpha,\beta)}(r)
%>   = \sqrt{n(n+\alpha+\beta+1)}\; \tilde{P}_{n-1}^{(\alpha+1,\beta+1)}(r). \f]
%>
%> @param[in] r     evaluation points within \f$[-1,1]\f$ (any shape)
%> @param[in] alpha Jacobi parameter, \f$\alpha > -1\f$
%> @param[in] beta  Jacobi parameter, \f$\beta > -1\f$
%> @param[in] n     polynomial order (non-negative integer)
%> @return dP derivative values, a numel(r)-by-1 column (same as the mex)
% ======================================================================
%> This function is part of the NDG-FEM software.
% ======================================================================

if ( n < 1 )
    dP = zeros( numel(r), 1 );
else
    dP = sqrt( n*(n + alpha + beta + 1.0) ) ...
        * JacobiP( r, alpha + 1.0, beta + 1.0, n - 1 );
end
end% func
