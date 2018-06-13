function gaugeValue = PostProcess_direction( obj )

Gauge_file =[ fileparts( mfilename('fullpath') ), '/tide/gaugepoint.txt'];
GaugePoint = load(Gauge_file);
xg = GaugePoint(:,2);
yg = GaugePoint(:,1);

pos = Analysis2d( obj, xg, yg );

Ng = numel(xg);
Nout = obj.outputFile.outputStep;
gaugeValue = zeros( Ng, Nout );

for i = 1:obj.outputFile.outputStep;
    temp(:,:) = pos.InterpOutputResult2GaugePoint( i );
    
    U = temp(:,2)./temp(:,1);
    V = temp(:,3)./temp(:,1);
    
    for j = 1 : pos.Ng
    
    if ( U(j)<0 && V(j)>=0)
        gaugeValue(j, i) = ...
            450 - atan2d(V(j), U(j));
    else
        gaugeValue(j, i) = ...
            90 - atan2d(V(j), U(j));
    end
    
    end
       
end


end

