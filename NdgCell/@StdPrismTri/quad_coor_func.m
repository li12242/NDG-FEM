function [ Nq, rq, sq, tq, wq ] = quad_coor_func( obj, N, N2 )
%> @brief Quadrature nodes and weights on the triangular prism,
%> as the tensor product of the collapsed-square triangle rule
%> (ndgcell.triquad) and the 1D LGL rule.

% horizontal triangle quadrature
qOrd = N+1;
[X,Y,Wx,Wy] = ndgcell.triquad(qOrd, [-1 -1; 1 -1; -1 1]);
rq1 = X(:);
sq1 = Y(:);
wq1 = Wx * Wy.';
wq1 = wq1(:);
Nq1 = numel(rq1);

% vertical 1D LGL quadrature
Nq2 = N2 + 1;
[ tq2, wq2 ] = zwglj( Nq2 );

Nq = Nq1 * Nq2;
r = repmat(rq1, 1, Nq2);
s = repmat(sq1, 1, Nq2);
wq1 = repmat(wq1, 1, Nq2);
t = repmat(tq2', Nq1, 1 );
wq2 = repmat(wq2', Nq1, 1 );
wq = wq1 .* wq2;

rq = r(:); sq = s(:); tq = t(:); wq = wq(:);
end
