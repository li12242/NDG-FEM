function [ field3d ] = Extend2dField( obj, field2d )
%EXTEND2DFIELD Summary of this function goes here
%   Detailed explanation goes here

K2d = obj.mesh2d.K;
field3d = zeros( obj.cell.Np, obj.K );

sk = ( (1:K2d) - 1 ) * obj.Nz + 1;
for n = 1 : obj.Nz
    field3d( :, sk ) = repmat( field2d, obj.cell.Npz, 1 );
    sk = sk + 1;
end

end
