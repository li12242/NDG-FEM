function [u] = Advec1D(u, FinalTime)

% function [u] = Advec1D(u, FinalTime)
% Purpose  : Integrate 1D advection until FinalTime starting with
%            initial the condition, u

% Globals1D;
time = 0;

% Runge-Kutta residual storage 
resu = scalar_type('resu', u.mesh);
resu = resu.allocate(size(u.val));

% compute time step size
xmin = min(abs(u.mesh.x(1,:)-u.mesh.x(2,:)));
CFL=0.75; dt   = CFL/(2*pi)*xmin; dt = .5*dt;
Nsteps = ceil(FinalTime/dt); dt = FinalTime/Nsteps; 

% advection speed
a = 2*pi;

% outer time step loop
for tstep=1:Nsteps
    for INTRK = 1:5
        timelocal = time + u.mesh.rk4c(INTRK)*dt;
        [rhsu] = AdvecRHS1D(u, timelocal, a);
        resu.val = u.mesh.rk4a(INTRK)*resu.val + dt*rhsu;
        u.val = u.val + u.mesh.rk4b(INTRK)* resu.val;
    end;
    % Increment time
    time = time+dt;
end;
return
