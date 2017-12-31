function [ gaugeValue ] = interpolateOutputStepResultToGaugePoint( obj, xg, yg, zg, outputStep )
Ng = numel( xg );
gaugeValue = zeros(Ng, obj.Nvar);
fresult = obj.accessOutputResultAtStepNum( outputStep ); %
for m = 1:obj.Nmesh
    [ cellId, Vg ] = obj.meshUnion.accessGaugePointLocation( xg, yg, zg );
    for n = 1:Ng
        if (cellId(n) == 0)
            % the gauge point locates out of the mesh domain
            continue;
        end
        for fld = 1:obj.Nvar
            gaugeValue(n, fld) = Vg(n, :) * fresult{m}(:, cellId(n), fld);
        end
    end
end
end