function gaugeValue = PostProcess_speed( obj )

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

    gaugeValue(:,i) = sqrt(( temp(:,2) ./ temp(:,1)).^2 + ...
            ( temp(:,3) ./ temp(:,1) ).^2);
        
end

end

