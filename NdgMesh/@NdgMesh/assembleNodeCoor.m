function [x, y, z] = assembleNodeCoor( obj, vx, vy, vz )
x = obj.proj_vert2node(vx);
y = obj.proj_vert2node(vy);
z = obj.proj_vert2node(vz);
end% func