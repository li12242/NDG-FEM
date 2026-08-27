function P = JacobiP( r, alpha, beta, n )
%> @brief Evaluate the orthonormal Jacobi polynomial \f$ \tilde{P}_n^{(\alpha,\beta)}(r) \f$.
%>
%> Replaces the Polylib mex gateway (thirdParty/Polylib/JacobiP.c + polylib.c).
%> The polynomial is first evaluated with the standard three-term recurrence
%> (Abramowitz & Stegun) and then scaled by the orthonormalization factor
%> \f[ \sqrt{ \frac{ (2n+\alpha+\beta+1)\, \Gamma(n+\alpha+\beta+1)\, \Gamma(n+1) }
%>   { 2^{\alpha+\beta+1}\, \Gamma(n+\alpha+1)\, \Gamma(n+\beta+1) } }, \f]
%> such that \f$ \int_{-1}^{1} \tilde{P}_n^2 w^{(\alpha,\beta)} dr = 1 \f$
%> (for \f$\alpha=\beta=0\f$, \f$\tilde{P}_0 = 1/\sqrt{2}\f$).
%>
%> @param[in] r     evaluation points within \f$[-1,1]\f$ (any shape)
%> @param[in] alpha Jacobi parameter, \f$\alpha > -1\f$
%> @param[in] beta  Jacobi parameter, \f$\beta > -1\f$
%> @param[in] n     polynomial order (non-negative integer)
%> @return P polynomial values, a numel(r)-by-1 column (same as the mex)
%> @note This port uses the floating-point gamma() for the normalization and
%> therefore stays exact for n >= 13, where the original C code overflows an
%> int32 factorial (thirdParty/Polylib/polylib.c:904).
% ======================================================================
%> This function is part of the NDG-FEM software.
% ======================================================================

if ( n < 0 ) || ( n ~= fix(n) )
    error('JacobiP:invalidOrder', 'Order n must be a non-negative integer.');
end
if ( alpha <= -1 ) || ( beta <= -1 )
    error('JacobiP:invalidPar', 'alpha and beta must be greater than -1.');
end

r = r(:); % mex always returns a numel(r)-by-1 column
Np = numel(r);

% three-term recurrence for the unnormalized polynomial P_n^{(alpha,beta)}
apo = ones(Np, 1); % P_0
if n == 0
    P = apo;
else
    aps = 0.5*( alpha - beta + (alpha + beta + 2.0).*r ); % P_1
    if n == 1
        P = aps;
    else
        apb = alpha + beta;
        for k = 2:n
            a1 = 2.0*k*(k + apb)*(2.0*k + apb - 2.0);
            a2 = (2.0*k + apb - 1.0)*(alpha*alpha - beta*beta)/a1;
            a3 = (2.0*k + apb - 2.0)*(2.0*k + apb - 1.0)*(2.0*k + apb)/a1;
            a4 = 2.0*(k + alpha - 1.0)*(k + beta - 1.0)*(2.0*k + apb)/a1;
            P = (a2 + a3.*r).*aps - a4*apo;
            apo = aps;
            aps = P;
        end
    end
end

% orthonormalization factor (gamma() instead of the C int32 factorial)
fac = (2.0*n + alpha + beta + 1.0)/2^(alpha + beta + 1) ...
    * gamma(n + alpha + beta + 1)*gamma(n + 1) ...
    / gamma(n + alpha + 1)/gamma(n + beta + 1);
P = sqrt(fac)*P;
end% func
