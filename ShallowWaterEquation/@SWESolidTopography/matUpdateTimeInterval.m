function dt = matUpdateTimeInterval( obj, fphys )

dt = nan;
for m = 1:obj.Nmesh
    dx = sqrt( obj.meshUnion(m).LAV );
    N = obj.meshUnion(m).cell.N;
    dtm = mxUpdateTimeInterval2d( obj.gra, N, dx, obj.meshUnion(m).EToR, fphys{m} );
    if ( dtm > 0 )
        dt = min(dt, dtm);
    end
end

end% func