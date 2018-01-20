function matEvaluateRK22( obj )

Nmesh = obj.Nmesh;
[rk4a, rk4b, rk4c] = GetRKParamter();

time = obj.getOption('startTime');
ftime = obj.getOption('finalTime');
resQ = cell( obj.Nmesh, 1 );
for n = 1:obj.Nmesh
    resQ{n} = zeros( obj.meshUnion(n).cell.Np, obj.meshUnion(n).K, obj.Nvar );
end
fphys = obj.fphys;
% init limiter and output file
hwait = waitbar(0,'Runing MatSolver....');
while( time < ftime )
    dt = obj.matUpdateTimeInterval( fphys )*0.5;
    if( time + dt > ftime )
        dt = ftime - time;
    end
    
    for intRK = 1:3
        tloc = time + rk4c(intRK) * dt;
        obj.matUpdateExternalField( tloc, fphys );
        obj.matEvaluateRHS( fphys );
        
        for n = 1:Nmesh
            resQ{n} = rk4a(intRK)*resQ{n} + dt*obj.frhs{n};
            fphys{n}(:,:, obj.varFieldIndex) ...
                = fphys{n}(:,:, obj.varFieldIndex) + rk4b(intRK)*resQ{n};
        end
        
        fphys = obj.matEvaluateLimiter( fphys );
        fphys = obj.matEvaluatePostFunc( fphys );
        
    end
%     for m = 1:obj.Nmesh
%         obj.meshUnion(m).draw( fphys{m}(:,:,1) );
%     end
%     drawnow;
    fprintf('processing %f ...\n', time/ftime);
    time = time + dt;
    obj.matUpdateOutputResult( time, fphys );
    %waitbar( time/ftime, hwait, 'Runing MatSolver....');
end
hwait.delete();
obj.matUpdateFinalResult( time, fphys );
obj.fphys = fphys;
end

function [rk4a, rk4b, rk4c] = GetRKParamter()
rk4a = [ 0.0, 0.0, -1 ];
rk4b = [ 0.5, 0.5, 1/3 ];
rk4c = [ 0.0, 0.5, 1 ];
end

