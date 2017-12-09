function [ f_Q ] = extval( obj, time )
%EXTVAL Summary of this function goes here
%   Detailed explanation goes here

w  = sqrt(8*obj.gra.*obj.a);
temp = obj.X+obj.Y*cos(w*time);

r2 = obj.mesh.x.^2 + obj.mesh.y.^2;
h = 1./temp + obj.a*(obj.Y^2 - obj.X^2)*r2./temp.^2;
h(h<0) = 0; f_Q(:,:,1) = h;

u = - obj.Y*w*sin(w*time)./temp.*obj.mesh.x./2;
v = - obj.Y*w*sin(w*time)./temp.*obj.mesh.y./2;

f_Q(:, :, 2) = u.*h;
f_Q(:, :, 3) = v.*h;

end

