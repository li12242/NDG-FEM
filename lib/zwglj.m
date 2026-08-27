function [ z, w ] = zwglj( np )
%> @brief Gauss-Lobatto-Legendre (GLL) points and weights on \f$[-1,1]\f$.
%>
%> Replaces the Polylib mex gateway (thirdParty/Polylib/zwglj.c) which
%> hard-codes \f$\alpha=\beta=0\f$. The np points always include the two
%> endpoints \f$\pm 1\f$ in ascending order; the np==1 degenerate case
%> returns \f$ z=0, w=2 \f$ (same as the mex).
%>
%> The interior np-2 points are the zeros of \f$P_{np-2}^{(1,1)}(z)\f$,
%> computed by the Golub-Welsch method (eigenvalues of the symmetric
%> tridiagonal Jacobi matrix, mirroring polylib.c:JacZeros/TriQL).
%> The weights follow the classical LGL formula
%> \f[ w_i = \frac{2}{np\,(np-1)}\, P_{np-1}(z_i)^{-2}, \f]
%> with \f$P_{np-1}\f$ the unnormalized Legendre polynomial.
%>
%> @param[in] np number of points (>= 1)
%> @return z point coordinates, np-by-1 column, ascending
%> @return w quadrature weights, np-by-1 column
%> @note Exactly two output arguments are required, same as the mex.
% ======================================================================
%> This function is part of the NDG-FEM software.
% ======================================================================

if ( nargout ~= 2 )
    error('zwglj:nargout', ...
        '[z,w] = zwglj(np) requires exactly two output arguments.');
end
if ( np < 1 ) || ( np ~= fix(np) )
    error('zwglj:invalidNp', 'np must be a positive integer.');
end

if ( np == 1 )
    z = 0; w = 2.0;
    return;
end

% interior nodes: zeros of P_{np-2}^{(1,1)} via Golub-Welsch
n = np - 2;        % number of interior points
alpha = 1.0; beta = 1.0; % raised parameters (alpha', beta' of polylib)
apb = alpha + beta;
if n >= 1
    a = zeros(n, 1); % diagonal entries, vanish for alpha == beta
    j = (1:(n-1)).'; % off-diagonal index, j = i+1 in polylib:JacZeros
    apbi = 2.0*j + apb;
    b = sqrt( 4.0*j.*(j + alpha).*(j + beta).*(j + apb) ...
        ./ ( (apbi.^2 - 1.0).*apbi.^2 ) );
    J = diag(a) + diag(b, 1) + diag(b, -1);
    zi = eig(J);
    zi = sort(zi); % ascending, mirroring TriQL's final selection sort
    z = [-1; zi; 1];
else
    z = [-1; 1];
end

% weights: fac/P_{np-1}(z)^2, fac = 2/(np*(np-1)) for alpha = beta = 0;
% the general endpoint factors (beta+1)/(alpha+1) are 1 here (polylib.c:288-297)
N = np - 1;
if N == 0
    P = ones(np, 1); %#ok<NASGU> % defensive, unreachable (np >= 2)
else
    Pm1 = ones(np, 1);
    P = z;
    for k = 2:N
        Pn = ( (2.0*k - 1.0)*z.*P - (k - 1.0)*Pm1 )/k;
        Pm1 = P;
        P = Pn;
    end
end
w = ( 2.0/( np*(np - 1.0) ) )./( P.*P );
end% func
