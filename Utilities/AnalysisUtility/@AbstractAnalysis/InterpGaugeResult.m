function data = InterpGaugeResult( obj, fphys )

Nfield = size( fphys{1}, 3 );
data = zeros( obj.Ng, Nfield );
for n = 1 : obj.Ng
    meshId = obj.gaugeMesh(n);
    cellId = obj.gaugeCell(n);
    Vg = obj.Vg{n};
    
    for fld = 1:Nfield
        data(n, fld) = Vg * fphys{meshId}(:, cellId, fld);
    end
end
end