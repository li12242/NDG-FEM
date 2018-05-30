function [ fdr, fds ] = derivative_orthogonal_func(Nh, ind, r, s)

    % transform the index to two indexes.
    i = mod( ind - 1, Nh + 1 );
    j = floor( (ind - 1) / (Nh + 1) );

    % calculate the derivative basis function values.
    fdr = GradJacobiP( r(:), 0, 0, i ).*JacobiP( s(:), 0, 0, j );
    fds = JacobiP( r(:), 0, 0, i ).*GradJacobiP( s(:), 0, 0, j );
end