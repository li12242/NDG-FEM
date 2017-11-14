function [ gaugeValue ] = interpolateOutputResultToGaugePoint( obj, xg, yg, zg )
Ntime = obj.Nt;
Ng = numel( xg );
gaugeValue = zeros(Ng, obj.Nfield, Ntime);
for m = 1:obj.Nmesh
    [ cellId, Vg ] = obj.meshUnion.accessGaugePointLocation( xg, yg, zg );
    % get the output file name
    Np = obj.meshUnion(m).cell.Np;
    K = obj.meshUnion(m).K;
    filename = [obj.getOption('outputNetcdfCaseName'), '.', ...
        num2str(m), '-', num2str(obj.Nmesh), '.nc'];
    % loop over all the output step
    for t = 1:Ntime
        fresult = ncread(filename, 'fphys', ...
            [1, 1, 1, t], [Np, K, obj.Nfield, 1]);
        for n = 1:Ng
            if (cellId(n) == 0)
                % the gauge point locates out of the mesh domain
                continue;
            end
            for fld = 1:obj.Nfield
                gaugeValue(n, fld, t) = Vg(n, :) * fresult(:, cellId(n), fld);
            end
        end
    end
end
end% func