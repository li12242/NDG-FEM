function matEvaluateEuler( obj )

Nmesh = obj.Nmesh;

time = obj.getOption('startTime');
ftime = obj.getOption('finalTime');

fphys = obj.fphys;
% init limiter and output file
hwait = waitbar(0,'Runing MatSolver....');
while( time < ftime )
    dt = obj.matUpdateTimeInterval( fphys )*0.5;
    if( time + dt > ftime )
        dt = ftime - time;
    end
    
    tloc = time + dt;
    obj.matUpdateExternalField( tloc, fphys );
    obj.matEvaluateRHS( fphys );
    
    for n = 1:Nmesh
        fphys{n}(:,:, obj.varFieldIndex) ...
            = fphys{n}(:,:, obj.varFieldIndex) + dt*obj.frhs{n};
    end
    
    fphys = obj.matEvaluateLimiter( fphys );
    fphys = obj.matEvaluatePostFunc( fphys );
    
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

