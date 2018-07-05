function matUpdateExternalField( obj, time, fphys2d, fphys3d )

delta = obj.externTimeInterval;

% find step
s1 = floor( time/delta ) + 1;
s2 = s1 + 1;
alpha1 = ( delta * (s2 - 1) - time ) / delta;
alpha2 = ( time - delta * (s1 - 1) ) / delta;

eta = obj.etaExt(s1) * alpha1 + obj.etaExt(s2) * alpha2;
for m = 1:obj.Nmesh
    obj.fext2d{m}(:, :, 1) = eta;
%     obj.fext2d{m}(:, :, 2) = ;

%     edge3d = obj.mesh3d(m).BoundaryEdge;
    obj.fext3d{m}(:, :, 7) = eta;
%     obj.fext3d{m}(:, :, 1) = - obj.a * sgH / obj.H / ...
%         cos( obj.omega * obj.ChLength / sgH ) .* ...
%         sin( obj.omega * ( obj.ChLength - edge3d.xb ) / sgH ) * ...
%         sin( obj.omega * time );
    
end

end
