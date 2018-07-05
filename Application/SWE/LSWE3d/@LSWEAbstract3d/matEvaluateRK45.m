function matEvaluateRK45( obj )

Nmesh = obj.Nmesh;
[rk4a, rk4b, rk4c] = GetRKParamter();

time = obj.startTime;
ftime = obj.finalTime;

resQ2d = cell( obj.Nmesh, 1 );
resQ3d = cell( obj.Nmesh, 1 );
for n = 1:obj.Nmesh
    Np2 = obj.mesh2d(n).cell.Np;
    K2 = obj.mesh2d(n).K;
    resQ2d{n} = zeros( Np2, K2, 1 );
    
    Np3 = obj.mesh3d(n).cell.Np;
    K3 = obj.mesh3d(n).K;
    resQ3d{n} = zeros( Np3, K3, 2 );
end

fphys2d = obj.fphys2d;
fphys3d = obj.fphys3d;

visual = Visual2d( obj.mesh2d );

dt = obj.dt;
hwait = waitbar(0,'Runing MatSolver....');
% try
while( time < ftime )
    %     dt = obj.matUpdateTimeInterval( fphys2d );
    if( time + dt > ftime )
        dt = ftime - time;
    end
    
    for intRK = 1:5
        tloc = time + rk4c( intRK ) * dt;
        obj.matUpdateExternalField( tloc, fphys2d, fphys3d );
        obj.matEvaluateRHS( fphys2d, fphys3d );
        
        for n = 1:Nmesh
            resQ2d{n} = rk4a( intRK ) * resQ2d{n} + dt * obj.frhs2d{n};
            resQ3d{n} = rk4a( intRK ) * resQ3d{n} + dt * obj.frhs3d{n};
            
            fphys2d{n}(:,:,1) = fphys2d{n}(:,:,1) + rk4b(intRK) * resQ2d{n};
            fphys3d{n}(:,:,1:2) = fphys3d{n}(:,:,1:2) + rk4b(intRK) * resQ3d{n};
            fphys3d{n}(:,:,3) = obj.matEvaluateVerticalVelocity( ...
                obj.mesh3d(n), fphys2d{n}, fphys3d{n} );
        end
        % fphys2d = obj.matEvaluateLimiter( fphys2d );
        % fphys2d = obj.matEvaluatePostFunc( fphys2d );
        % visual.drawResult( fphys2d{1}(:,:,1) );
        % figure; obj.mesh3d.drawHorizonSlice( fphys3d{1}(:, :, 1) )
    end
    
    % visual.drawResult( fphys2d{1}(:,:,1) );
    % obj.drawVerticalSlice( 20, 1, fphys3d{1}(:, :, 3) * 1e7 );
    time = time + dt;
    obj.matUpdateOutputResult( time, fphys2d, fphys3d );
    
    timeRatio = time / ftime;
    waitbar( timeRatio, hwait, ['Runing MatSolver ', num2str( timeRatio ), '....']);
end
hwait.delete();
obj.fphys2d = fphys2d;
obj.fphys3d = fphys3d;

obj.outputFile.closeOutputFile();
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