function matEvaluateRK45( obj )

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
    dt = obj.matUpdateTimeInterval( fphys );
    if( time + dt > ftime )
        dt = ftime - time;
    end
    
    for intRK = 1:5
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
%     obj.draw( fphys );
%     obj.meshUnion(1).draw( fphys{1}(:,:,1) + fphys{1}(:,:,4) );
%     obj.meshUnion(1).draw( fphys{1}(:,:,1) );
%     obj.meshUnion(1).draw( fphys{1}(:,:,1) );
%     drawnow;
    
    time = time + dt;
    obj.matUpdateOutputResult( time, fphys );
    %fprintf('processing %f ...\n', time/ftime);
    timeRatio = time/ftime;
    waitbar( timeRatio, hwait, ...
        ['Runing MatSolver ', num2str(timeRatio), '....']);
end
hwait.delete();
obj.matUpdateFinalResult( time, fphys );
obj.fphys = fphys;
end

function [rk4a, rk4b, rk4c] = GetRKParamter()
rk4a = [            0.0 ...
    -567301805773.0/1357537059087.0 ...
    -2404267990393.0/2016746695238.0 ...
    -3550918686646.0/2091501179385.0  ...
    -1275806237668.0/842570457699.0];
rk4b = [ 1432997174477.0/9575080441755.0 ...
    5161836677717.0/13612068292357.0 ...
    1720146321549.0/2090206949498.0  ...
    3134564353537.0/4481467310338.0  ...
    2277821191437.0/14882151754819.0];
rk4c = [             0.0  ...
    1432997174477.0/9575080441755.0 ...
    2526269341429.0/6820363962896.0 ...
    2006345519317.0/3224310063776.0 ...
    2802321613138.0/2924317926251.0];
end