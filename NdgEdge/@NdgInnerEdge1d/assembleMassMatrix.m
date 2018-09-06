function obj = assembleMassMatrix( obj )
%ASSEMBLEMASSMATRIX Summary of this function goes here
%   Detailed explanation goes here

cell = StdPoint( obj.mesh.cell.N );
obj.M = cell.M;
obj.Nfp = cell.Np;
obj.r = cell.r;

end

