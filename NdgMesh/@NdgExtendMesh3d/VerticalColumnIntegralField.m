function [ field2d ] = VerticalColumnIntegralField( obj, field3d )
%VERTICALINTEGRALFIELD Summary of this function goes here
%   Detailed explanation goes here

Np2 = obj.mesh2d.cell.Np;
K2d = obj.mesh2d.K;

fmod = obj.cell.V \ (obj.Jz .* field3d);
field2d = zeros( Np2, K2d );

sk = ( (1:K2d) - 1 ) * obj.Nz + 1;
for n = 1:obj.Nz
    field2d = field2d + fmod( 1:Np2, sk );
    sk = sk + 1;
end

% the mon coefficient need to devide the \varphi_{z,0}
field2d = obj.mesh2d.cell.V * ( field2d ./ obj.mesh2d.cell.V(:, 1) );
end
