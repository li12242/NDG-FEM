function c = TransDiffSolver(mesh, c, FinalTime)
% 1D Transient Diffusion Solver
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
resQ = zeros(size(c)); 
% compute time step size
xmin = min(abs(mesh.x(1,:)-mesh.x(2,:)));
% time step
CFL = 0.125; dt = (xmin)^2*CFL;

% intermediate variable
% q = zeros(c);

% outer time step loop 
while(time<FinalTime)
    for INTRK = 1:5
        rhsC = TransDiffRHS(mesh, c, time);%, invM);
        resQ = rk4a(INTRK)*resQ + dt*rhsC;
        c = c + rk4b(INTRK)*resQ;
%         c = Utilities.SlopeLimit1(mesh, c);
%         TranDiffPostProcess(mesh, c); drawnow;
    end% for
    % Increment time
    time = time+dt;
    fprintf(['Process: ', num2str(time/FinalTime*100), '%% \n'])
end% while
end% func