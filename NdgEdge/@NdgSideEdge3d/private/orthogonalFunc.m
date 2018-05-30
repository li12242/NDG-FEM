function f = orthogonalFunc( ind, Nh, r, s )
    i = mod( ind - 1, Nh + 1 );
    j = floor( (ind - 1) / (Nh + 1) );

    f = JacobiP( r(:), 0, 0, i ) .* JacobiP( s(:), 0, 0, j );
end