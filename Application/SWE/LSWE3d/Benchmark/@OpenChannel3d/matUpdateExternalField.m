function matUpdateExternalField( obj, time, fphys2d, fphys3d )
%MATUPDATEEXTERNALFIELD Summary of this function goes here
%   Detailed explanation goes here

eta = obj.a * cos( obj.omega * time );
sgH = sqrt( obj.gra * obj.H );

for m = 1:obj.Nmesh
    edge2d = obj.mesh2d(m).BoundaryEdge;
    obj.fext2d{m}(:, :, 1) = eta;
    obj.fext2d{m}(:, :, 2) = - obj.a * sgH / ...
        cos( obj.omega * obj.ChLength / sgH ) .* ...
        sin( obj.omega * ( obj.ChLength - edge2d.xb ) / sgH ) * ...
        sin( obj.omega * time );

    edge3d = obj.mesh3d(m).BoundaryEdge;
    obj.fext3d{m}(:, :, 7) = eta;
    obj.fext3d{m}(:, :, 1) = - obj.a * sgH / obj.H / ...
        cos( obj.omega * obj.ChLength / sgH ) .* ...
        sin( obj.omega * ( obj.ChLength - edge3d.xb ) / sgH ) * ...
        sin( obj.omega * time );
    
end
end
