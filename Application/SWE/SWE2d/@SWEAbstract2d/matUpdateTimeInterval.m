%> \brief determine the time interval
function [ dt ] = matUpdateTimeInterval( obj, fphys )
dt = nan;
for m = 1:obj.Nmesh
    dx = sqrt( obj.meshUnion(m).LAV );
    N = obj.meshUnion(m).cell.N;
    dtm = mxUpdateTimeInterval2d( obj.hmin, ...
        obj.gra, N, dx, obj.meshUnion(m).status, fphys{m} );
    
    if ( dtm > 0 )
        dt = min(dt, dtm * obj.cfl);
    end
end

end

