function gaugeValue = PostProcess_direction( obj, ncfilename )
sanyapos = NdgPostProcess(obj.meshUnion,ncfilename);
for t = 1:sanyapos.Nt
    field = sanyapos.accessOutputResultAtStepNum( t );
    
    for m = 1:sanyapos.Nmesh
        Uspeed{m}(:,:,t) = ( field{m}(:,:,2) ./ field{m}(:,:,1) );
        Vspeed{m}(:,:,t) = ( field{m}(:,:,3) ./ field{m}(:,:,1) );
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
        fresult1 = Uspeed{m}(:,:,t);
        fresult2 = Vspeed{m}(:,:,t);
        for n = 1:Ng
            if (cellId(n) == 0)
                % the gauge point locates out of the mesh domain
                continue;
            end
            UVValue{m}(n, t,1) = Vg(n, :) * fresult1(:, cellId(n));
            UVValue{m}(n, t,2) = Vg(n, :) * fresult2(:, cellId(n));
            if ( UVValue{m}(n, t,1)<0 && UVValue{m}(n, t,2)>=0)
                gaugeValue(n, t) = ...
                    450 - atan2d(UVValue{m}(n, t,2), UVValue{m}(n, t,1));
            else
                gaugeValue(n, t) = ...
                    90 - atan2d(UVValue{m}(n, t,2), UVValue{m}(n, t,1));
            end
            % gaugeValue(n, t,1) = Vg(n, :) * fresult1(:, cellId(n));
            % gaugeValue(n, t,2) =  Vg(n, :) * fresult2(:, cellId(n));
        end
    end
end

end

