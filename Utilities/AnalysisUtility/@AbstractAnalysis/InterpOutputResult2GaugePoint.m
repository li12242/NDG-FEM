function [ data ] = InterpOutputResult2GaugePoint( obj, step )
%INTERPOUTPUTRESULT2GAUGEPOINT Summary of this function goes here
%   Detailed explanation goes here

Nstep = obj.phys.outputFile(1).outputStep;
if( step > Nstep )
    error( 'Output step num is out of range.' );
end

data = zeros( obj.Ng, obj.phys.Nvar );
for m = 1 : obj.Nmesh
    if( obj.NgaugePerMesh(m) > 0 )
        gid = obj.gaugeIdPerMesh( 1:obj.gaugeIdPerMesh(obj.NgaugePerMesh(m)) );
        field = obj.phys.outputFile( m ).readOutputResult( step );
        
        for n = 1:obj.NgaugePerMesh(m)
            id = gid(n);
            cellId = obj.gaugeCell(id);
            Vg = obj.Vg{id};
            
            for fld = 1:obj.phys.Nvar
                data(id, fld) = Vg * field(:, cellId, fld);
            end
        end
    end
end


% for n = 1 : obj.Ng
%     meshId = obj.gaugeMesh(n);
%     cellId = obj.gaugeCell(n);
%     Vg = obj.Vg{n};
%     
%     Nstep = obj.phys.outputFile(meshId).outputStep;
%     if( step > Nstep )
%         error( 'Output step num is out of range.' );
%     end
%     
%     field = obj.phys.outputFile( meshId ).readOutputResult( step );
%     for fld = 1:obj.phys.Nvar
%         data(n, fld) = Vg * field(:, cellId, fld);
%     end
%     
% end

end
