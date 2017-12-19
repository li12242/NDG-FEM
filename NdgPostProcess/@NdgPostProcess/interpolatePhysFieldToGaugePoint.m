function [ gaugeValue ] = interpolatePhysFieldToGaugePoint( obj, fphys, xg, yg, zg )
Ng = numel( xg );
gaugeValue = zeros(Ng, obj.Nvar);

for m = 1:obj.Nmesh
    [ cellId, Vg ] = obj.meshUnion(m).accessGaugePointLocation( xg, yg, zg );
    for n = 1:Ng
        if (cellId(n) == 0)
            % the gauge point locates out of the mesh domain
            continue;
        end
        
        for fld = 1:obj.Nvar
            gaugeValue(n, fld) = Vg(n, :) * fphys{m}(:, cellId(n), fld);
        end
    end
end
end