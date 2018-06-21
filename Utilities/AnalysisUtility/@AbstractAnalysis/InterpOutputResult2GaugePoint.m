function [ data ] = InterpOutputResult2GaugePoint( obj, step )
%INTERPOUTPUTRESULT2GAUGEPOINT Summary of this function goes here
%   Detailed explanation goes here

data = zeros( obj.Ng, obj.phys.Nvar );
for n = 1 : obj.Ng
    meshId = obj.gaugeMesh(n);
    cellId = obj.gaugeCell(n);
    Vg = obj.Vg{n};
    
    Nstep = obj.phys.outputFile(meshId).outputStep;
    if( step > Nstep )
        error( 'Output step num is out of range.' );
    end
    
    field = obj.phys.outputFile( meshId ).readOutputResult( step );
    for fld = 1:obj.phys.Nvar
        data(n, fld) = Vg * field(:, cellId, fld);
    end
    
end

end
