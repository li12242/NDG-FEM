function [Q] = SWE2D(mesh,Q,FinalTime,fluxFunc)

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
     
time = 0;
% Runge-Kutta residual storage  
resQ = zeros(size(Q));

% compute time step size
xmin = min(sqrt(abs(mesh.x(1,:)-mesh.x(2,:)).^2 + ...
    abs(mesh.y(1,:)-mesh.y(2,:)).^2));

Nfp = @(x)mesh.Shape.nFaceNode;
Dr = @(x)mesh.Shape.Dr; 
Ds = @(x)mesh.Shape.Ds;
Fmat = @(x)mesh.Shape.FaceMassMatrixSmall;
invM = @(x)mesh.Shape.invM;

dt=0.001; CFL=0.010; 
% outer time step loop 
while(time<FinalTime)
    for INTRK = 1:5
%         timelocal = time + rk4c(INTRK)*dt;
        [rhsQ, lamda] = fluxFunc(mesh, Q, Nfp, Dr, Ds, Fmat, invM);
        resQ = rk4a(INTRK)*resQ + dt*rhsQ;
        Q = Q + rk4b(INTRK)*resQ;
%         plot(mesh.x, Q(:,:,1), 'r'); drawnow;
    end;
%     Q = real(Q);
%     for n = 1:2
%         Q(:,:,n)  = Utilities.SlopeLimit1(mesh, Q(:,:,n));
%     end
    plot3(mesh.x, mesh.y, Q(:,:,1), 'ro'); drawnow;
    % Increment time
    time = time+dt;
    dt = CFL/max(lamda(:))*xmin;
end
end