% max order of polymomials
N = 1; nElement = 11; x1 = -1000; x2 = 1000; % Parabolic Bowl
[Nv, VX, ~, EToV] = Utilities.Mesh.MeshGen1D(x1, x2, nElement);
BC = [2,1; 3,Nv];

line = StdRegions.Line(N);
mesh = MultiRegions.RegionLineBC(line, EToV, VX, BC);

a = 600; h0 = 10;
VB = h0.*(VX.^2./a^2 - 1);

vb = VB(EToV');
bedElevation = 0.5*((1-mesh.Shape.r)*vb(1,:) + (mesh.Shape.r+1)*vb(2,:));

% find the wet/dry interface
eta0 = -2;
index = find(bedElevation(1,:).*bedElevation(end,:)<0);
id = index(1); 
x1 = mesh.x(1, id); x2 = mesh.x(end, id); 
y1 = bedElevation(1, id); y2 = bedElevation(end, id);
xw = x2 + (eta0 - y2)*(x1 - x2)./(y1 - y2);
yw = y2 + (y1 - y2)*(xw - x2)./(x1 - x2);
volume = .5*(x2 - xw)*(eta0 - y2);
etaw = 0.5*(y1 + y2) + volume./(x2 - x1);

eta = eta0*ones(size(mesh.x));
depth = (bedElevation(end,id+1) - 0.5)*ones(size(mesh.x));

h = eta - bedElevation;
h(h<0) = 0;
h = Utilities.Limiter.Limiter1D.MinmodLinear(mesh,h);

% draw figure
close all; figure('Position', [220   705   568   280], 'Color', 'w');
plot(mesh.x(:, [id:id+1]), bedElevation(:, [id:id+1]), 'k.-', 'Markersize', 20); hold on
plot(mesh.x(:, id+1), eta(:, id+1), 'b.-', mesh.x(:, id), [etaw; etaw], 'r.-', 'Markersize', 20)
plot([xw, x2], [yw, yw], 'b.--', [x1, xw], [y1, yw], 'b.--', 'Markersize', 20)
plot([x1, x1], [depth(1), y1 + 0.5], 'k--', [x2, x2], [depth(1), y1 + 0.5], 'k--')
x3 = mesh.x(2, id+1);
plot([x3, x3], [depth(1), y1 + 0.5], 'k--')

xLims = get(gca, 'XLim'); yLims = get(gca, 'YLim'); 
plot(mesh.x(:), depth(:), 'k.-', xw, depth(1), 'k.', 'Markersize', 20)
t(1) = text(x1+10, y1+0.1, '${\it B}_{i-1/2}$'); 
t(2) = text(x2+10, y2+0.1, '${\it B}_{i+1/2}$');
t(3) = text(x1+5, etaw-0.9, '$\eta_{i-1/2}^+$');
t(4) = text(x2+5, etaw+0.8, '$\eta_{i+1/2}^-$');
t(5) = text(x2+5, eta0 - 0.7, '$\eta_{i+1/2}^+$');
t(6) = text(xw-10, depth(1) - 0.9, '${\it x}_{\omega}^{\ast}$');
t(7) = text((x1 + x2)./2, y1 - 0.5, '$\Omega_i$');
t(8) = text((x3 + x2)./2, y1 - 0.5, '$\Omega_{i+1}$');
set(t, 'Interpreter','latex', 'FontSize', 15)
set(t(7:8), 'FontSize', 20)

xlim(xLims); ylim([depth(1)-1.5, yLims(2)]);
set(gca, 'DataAspectRatio', [20, 1, 1])
set(gca, 'XTick', [], 'YTick', [])

% reconstruction result
eta1 = h + bedElevation;
eta1(1, id) = bedElevation(1, id);
figure('Position', [220   705   568   280], 'Color', 'w');
plot(mesh.x(:, [id:id+1]), bedElevation(:, [id:id+1]), 'k.-', 'Markersize', 20); hold on
plot(mesh.x(:, id+1), eta(:, id+1), 'b.-', mesh.x(:, id), eta1(:, id), 'r.-', 'Markersize', 20)
plot([xw, x2], [yw, yw], 'b.--', [x1, xw], [y1, yw], 'b.--', 'Markersize', 20)
plot([x1, x1], [depth(1), y1 + 0.5], 'k--', [x2, x2], [depth(1), y1 + 0.5], 'k--')
x3 = mesh.x(2, id+1);
plot([x3, x3], [depth(1), y1 + 0.5], 'k--')


xLims = get(gca, 'XLim'); yLims = get(gca, 'YLim'); 
plot(mesh.x(:), depth(:), 'k.-', xw, depth(1), 'k.', 'Markersize', 20)
t(1) = text(x1+10, y1+0.1, '${\it B}_{i-1/2}$'); t(2) = text(x2+10, y2+0.1, '${\it B}_{i+1/2}$');
t(3) = text(x2+5, eta1(2, id)-0.2, '$\eta_{i+1/2}^-$');
t(4) = text(x2+5, eta0 + 0.9, '$\eta_{i+1/2}^+$');
t(5) = text(xw-10, depth(1) - 0.9, '${\it x}_{\omega}^{\ast}$');
t(7) = text((x1 + x2)./2, y1 - 0.5, '$\Omega_i$');
t(8) = text((x3 + x2)./2, y1 - 0.5, '$\Omega_{i+1}$');
set(t, 'Interpreter','latex', 'FontSize', 15)
xlim(xLims); ylim([depth(1)-1.5, yLims(2)]);
set(gca, 'DataAspectRatio', [20, 1, 1])
set(gca, 'XTick', [], 'YTick', [])
set(t(7:8), 'FontSize', 20)
