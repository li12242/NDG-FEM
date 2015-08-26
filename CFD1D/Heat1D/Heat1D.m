function [u,time] = Heat1D(u,FinalTime,Element, Mesh)

% function [u] = Heat1D(u,FinalTime)
% Purpose  : Integrate 1D heat equation until 
%            FinalTime starting with initial condition, u.

% Globals1D;
Np = Element.Np; K = Mesh.K;

time = 0;

% Runge-Kutta residual storage  
resu = zeros(Np, K); 

% compute time step size
xmin = min(abs(Mesh.x(1,:)-Mesh.x(2,:)));
CFL=0.25;dt   = CFL*(xmin)^2;
Nsteps = ceil(FinalTime/dt); dt = FinalTime/Nsteps;

% outer time step loop 
for tstep=1:Nsteps
  for INTRK = 1:5
    timelocal = time + Mesh.rk4c(INTRK)*dt;        

    % compute right hand side of 1D advection equations
    [rhsu] = HeatCRHS1D(u,timelocal,Element,Mesh);

    % initiate and increment Runge-Kutta residuals
    resu = Mesh.rk4a(INTRK)*resu + dt*rhsu;  
    
    % update fields
    u = u+Mesh.rk4b(INTRK)*resu;
  end;
  % Increment time
  time = time+dt;
end
return
