function refine_mesh( obj )
%REFINE_MESH Summary of this function goes here
%   Detailed explanation goes here

obj.mesh = obj.mesh.refine(1);
obj.init;

end

