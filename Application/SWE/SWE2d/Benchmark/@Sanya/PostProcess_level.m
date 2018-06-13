function gaugeValue = PostProcess_level( obj )

Gauge_file =[ fileparts( mfilename('fullpath') ), '/tide/gaugepoint.txt'];
GaugePoint = load(Gauge_file);
xg = GaugePoint(:,2);
yg = GaugePoint(:,1);

pos = Analysis2d( obj, xg, yg );
field = pos.InterpGaugeResult( obj.fphys );

Ng = numel(xg);
Nout = obj.outputFile.outputStep;
gaugeValue = zeros( Ng, Nout );
for i = 1:Nout
    temp = pos.InterpOutputResult2GaugePoint( i );
    gaugeValue(:,i) = temp(:,1) + field(:,4);
end

end





