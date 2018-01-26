function [ gaugeValue ] = interpolateOutputResultToGaugePoint( obj, xg, yg, zg )
Ntime = obj.Nt;
Ng = numel( xg );
gaugeValue = zeros(Ng, obj.Nvar, Ntime);
for m = 1:obj.Nmesh
    [ cellId, Vg ] = obj.meshUnion.accessGaugePointLocation( xg, yg, zg );
    % get the output file name
    Np = obj.meshUnion(m).cell.Np;
    K = obj.meshUnion(m).K;
    % loop over all the output step
    for t = 1:Ntime
        fresult = ncread(obj.outputFile{m}, 'fphys', ...
            [1, 1, 1, t], [Np, K, obj.Nvar, 1]);
        for n = 1:Ng
            if (cellId(n) == 0)
                % the gauge point locates out of the mesh domain
                continue;
            end
            for fld = 1:obj.Nvar
                gaugeValue(n, fld, t) = Vg(n, :) * fresult(:, cellId(n), fld);
            end
        end
    end
end
end% func