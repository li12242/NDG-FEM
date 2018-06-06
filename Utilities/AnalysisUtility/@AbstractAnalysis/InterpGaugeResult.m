function InterpGaugeResult( obj, fieldId )
    
    for n = 1 : obj.Ng
        figure(n)
        meshId = obj.gaugeMesh(n);
        cellId = obj.gaugeCell(n);
        Vg = obj.Vg{n};

        Nstep = obj.phys.outputFile(meshId).outputStep;
        time = obj.phys.outputFile(meshId).outputTime;
        data = zeros( Nstep, 1 );
        for i = 1 : Nstep
            field = obj.phys.outputFile( meshId ).readOutputResult( i );
            data(i) = Vg * field(:, cellId, fieldId);
        end
        
        plot( time, data ); 
    end
end