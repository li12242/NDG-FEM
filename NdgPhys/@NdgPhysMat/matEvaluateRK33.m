function matEvaluateRK33( obj )

Nmesh = obj.Nmesh;
[rk0, rk1, rk2, rkt] = GetRKParamter();

time = obj.getOption('startTime');
ftime = obj.getOption('finalTime');
fphys0 = cell( obj.Nmesh, 1 );
for n = 1:obj.Nmesh
    fphys0{n} = zeros( obj.meshUnion(n).cell.Np, obj.meshUnion(n).K, obj.Nvar );
end
fphys = obj.fphys;
% init limiter and output file
hwait = waitbar(0,'Runing MatSolver....');
while( time < ftime )
    dt = obj.matUpdateTimeInterval( fphys ) * 0.5;
    if( time + dt > ftime )
        dt = ftime - time;
    end
    
    for n = 1:obj.Nmesh
        fphys0{n} = fphys{n};
    end
    
    for intRK = 1:3
        tloc = time + rkt(intRK) * dt;
        obj.matUpdateExternalField( tloc, fphys );
        obj.matEvaluateRHS( fphys );
        
        for n = 1:Nmesh
            fphys{n}(:,:, obj.varFieldIndex) ...
                = rk0(intRK)*fphys0{n}(:,:, obj.varFieldIndex)...
                + rk1(intRK)*fphys{n}(:,:, obj.varFieldIndex) ...
                + rk2(intRK)*obj.frhs{n};
        end
        
        fphys = obj.matEvaluateLimiter( fphys );
        fphys = obj.matEvaluatePostFunc( fphys );
        
    end
    for m = 1:obj.Nmesh
        obj.meshUnion(m).draw( fphys{m}(:,:,1) );
    end
    drawnow;
    
    time = time + dt;
    obj.matUpdateOutputResult( time, fphys );
    waitbar( time/ftime, hwait, 'Runing MatSolver....');
end
hwait.delete();
obj.matUpdateFinalResult( time, fphys );
obj.fphys = fphys;
end

function [rk0, rk1, rk2, rkt] = GetRKParamter()
rk0 = [ 1.0, 3/4, 1/3 ];
rk1 = [ 0.0, 1/4, 2/3 ];
rk2 = [ 1.0, 1/4, 2/3 ];
rkt = [ 0.0, 1.0, 1.0];
end

