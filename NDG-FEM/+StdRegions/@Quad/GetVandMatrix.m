function Va = GetVandMatrix(obj, order, r, s)
% GETVANDMATRIX ¼ÆËãVandermonde¾ØÕó
% Vandermonde¾ØÕóÔªËØÂú×ã V_{ij} = \psi_j(\xi_i)

Va = zeros(numel(r), obj.nNode);

for ind = 1:obj.nNode
    Va(:, ind) = OrthogonalFun(r,s,order,ind);
end