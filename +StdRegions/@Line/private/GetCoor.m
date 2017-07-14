function r = GetCoor(nOrder)
% get Gauss Lobatto jacobi points at [-1, 1]
[r,~] = Polylib.zwglj(nOrder+1);
end