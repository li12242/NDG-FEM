function refine_mesh( obj, refine_level )
%REFINE_MESH Summary of this function goes here
%   Detailed explanation goes here

obj.mesh = obj.mesh.refine(refine_level);
obj.init;

end

