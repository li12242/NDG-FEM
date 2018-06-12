function data = InterpGaugeResult( obj, fphys )
    
    data = zeros( obj.Ng, obj.Nfield );
    for n = 1 : obj.Ng
        meshId = obj.gaugeMesh(n);
        cellId = obj.gaugeCell(n);
        Vg = obj.Vg{n};

        for fld = 1:obj.Nfield
            data(n, fld) = Vg * fphys{meshId}(:, cellId, fld);
        end
    end
end