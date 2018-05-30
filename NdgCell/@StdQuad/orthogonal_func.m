function [ fval ] = orthogonal_func(obj, N, ind, r, s, t)

i = mod( ind - 1, N + 1 );
j = floor( (ind - 1) / (N + 1) );
fval = JacobiP( r(:), 0, 0, i ).*JacobiP( s(:), 0, 0, j );

end

