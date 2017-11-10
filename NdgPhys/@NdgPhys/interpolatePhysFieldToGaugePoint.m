function [ gaugeValue ] = interpolatePhysFieldToGaugePoint( obj, xg, yg, zg )
Ng = numel( xg );
gaugeValue = zeros(Ng, obj.Nfield);

for m = 1:obj.Nmesh
    [ cellId, Vg ] = obj.meshUnion.accessGaugePointLocation( xg, yg, zg );
    for n = 1:Ng
        if (cellId(n) == 0)
            % the gauge point locates out of the mesh domain
            continue;
        end
        
        for fld = 1:obj.Nfield
            gaugeValue(n, fld) = Vg(n, :) * obj.fphys{m}(:, cellId(n), fld);
        end
    end
end
end