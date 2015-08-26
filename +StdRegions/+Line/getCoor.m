function r = getCoor(nOrder)
% get Gauss Lobatto jacobi points at [-1, 1]
[r,~] = Polylib.zwglj(nOrder+1); r=r';
end