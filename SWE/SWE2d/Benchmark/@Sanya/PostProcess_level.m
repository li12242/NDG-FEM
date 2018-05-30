function gaugeValue = PostProcess_level( obj, ncfilename )
sanyapos = NdgPostProcess(obj.meshUnion,ncfilename);
for t = 1:sanyapos.Nt
    field = sanyapos.accessOutputResultAtStepNum( t );
    
    for m = 1:sanyapos.Nmesh
        TideElevation{m}(:,:,t) = ( field{m}(:,:,1) + ...
            obj.fphys{1}(:,:,4) );
    end
end

Ntime = sanyapos.Nt;

Gauge_file =[pwd, '/SWE2d/@Sanya/tide/gaugepoint.txt'];
fp = fopen(Gauge_file,'r');
GaugePoint = fscanf(fp, '%f\n', [3, inf]);
fclose(fp); 
GaugePoint = GaugePoint';
xg = GaugePoint(:,2);
yg = GaugePoint(:,1);
zg = GaugePoint(:,3);

Ng = numel( xg );
gaugeValue = zeros(Ng, Ntime);

    for m = 1:sanyapos.Nmesh
        [ cellId, Vg ] = sanyapos.meshUnion.accessGaugePointLocation( xg, yg, zg );
        
        for t = 1:Ntime
            fresult = TideElevation{m}(:,:,t);
            for n = 1:Ng
                if (cellId(n) == 0)
                    % the gauge point locates out of the mesh domain
                    continue;
                end
                gaugeValue(n, t) = Vg(n, :) * fresult(:, cellId(n));
            end
        end
    end



end





