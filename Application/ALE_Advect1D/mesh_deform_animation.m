% deform amination
x1 = 0; x2 = 2*pi; % space domain
N = 3; nElement = 20; delta_x = (x2 - x1)/nElement;
[Nv, VX, K, EToV] = Utilities.Mesh.MeshGen1D(x1, x2, nElement);

FinalTime = 10; nT = 40;
t = linspace(0, FinalTime, nT);
w = 2*pi/(FinalTime/3); % period T = finalTime/3 = 3.3s
x = VX;
for it = 1:nT
    temp = VX + delta_x/4.*sin( w*t(it) + VX );
    x(2:end-1) = temp(2:end-1);
    plot(x, zeros(size(x)), 'ro');
    drawnow; pause(0.1);
end
