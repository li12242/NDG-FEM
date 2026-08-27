function [ Nq, rq, sq, tq, wq ] = quad_coor_func( obj, N )
%> @brief Quadrature nodes and weights on the reference triangle.
%> Uses the collapsed-square Gauss rule (ndgcell.triquad) of order N+1.
qOrd = N+1;
[X,Y,Wx,Wy] = ndgcell.triquad(qOrd, [-1 -1; 1 -1; -1 1]);
rq = X(:);
sq = Y(:);
wq = Wx * Wy.';
wq = wq(:);

Nq = numel(rq);
tq = zeros(Nq, 1);
end% func
