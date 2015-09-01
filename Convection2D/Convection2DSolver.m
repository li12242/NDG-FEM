function var = Convection2DSolver(mesh, var, FinalTime, Speed)
% 2D convection problem
% RK time stepping

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
     
time = 0;
resVar = zeros(size(var));
% compute time step size
xmin = min(sqrt((mesh.x(1,:)-mesh.x(2,:)).^2 + (mesh.y(1,:) - mesh.y(2,:)).^2 ));
CFL=0.30;  un = sqrt(sum(Speed.^2));
dt = CFL/un*xmin;

while(time < FinalTime)
    time = time+dt;
    fprintf('Processing: %f ...\n', time./FinalTime)
    for INTRK = 1:5
        rhsVar = Convection2DRHS(mesh, var, time, Speed);
        resVar = rk4a(INTRK)*resVar + dt*rhsVar;
        var = var + rk4b(INTRK)*resVar;
        plot3([mesh.x; mesh.x(1,:)], [mesh.y; mesh.y(1,:)], [var; var(1,:)], '.-'); drawnow
    end% for
    

end% while


end% func