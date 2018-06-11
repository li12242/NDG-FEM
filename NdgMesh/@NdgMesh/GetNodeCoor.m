function obj = GetNodeCoor( obj )

    obj.x = obj.proj_vert2node(obj.vx);
obj.y = obj.proj_vert2node(obj.vy);
obj.z = obj.proj_vert2node(obj.vz);

end% func