function [ fval ] = orthogonal_func(obj, N, ind, r, s, t)

[ i,j ] = trans_ind(N, ind);
fval = JacobiP(r(:), 0, 0, i).*JacobiP(s(:), 0, 0, j);

end

