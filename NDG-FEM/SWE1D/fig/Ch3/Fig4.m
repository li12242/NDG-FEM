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

% element refinement
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
depth = (eta0 - 3)*ones(size(mesh.x));

h = eta - bedElevation;
h(h<0) = 0;
h = Utilities.Limiter.Limiter1D.MinmodLinear(mesh,h);
eta1 = h + bedElevation;
eta1(1, id) = bedElevation(1, id);
% draw original
close all; figure('Position', [220   705   568   280], 'Color', 'w');
plot(mesh.x(:,id), bedElevation(:, id), 'k.-', 'MarkerSize', 20); hold on;
set(gca, 'XTick', [], 'YTick', []); yLims = get(gca, 'YLim');
plot(mesh.x(:, id), depth(:, id), 'k.-', xw, depth(1), '.', 'MarkerSize', 20);
plot(mesh.x(:, id), eta1(:, id), 'b.-', 'Markersize', 20)
plot([x1, x1], [depth(1), y1 + 0.5], 'k--', [x2, x2], [depth(1), y1 + 0.5], 'k--')

ylim([depth(1)-1, yLims(2)])
t(1) = text(xw-5, depth(1) - 0.5, '${\it x}_{\omega}^{\ast}$');
t(2) = text(x1, depth(1) - 0.5, '${\it x}_{i-1/2}$');
t(3) = text(x2-7, depth(1) - 0.5, '${\it x}_{i+1/2}$');
t(4) = text(x1+3, y1-0.9, '$B_{i-1/2}$');
t(5) = text(x2-20, y2-0.3, '$B_{i+1/2}$');
t(6) = text((x1 + x2)./2, y1 - 0.5, '$\Omega_i$');
set(t, 'Interpreter','latex', 'FontSize', 15)
set(t(6), 'FontSize', 20)

% refinement elements
figure('Position', [220   705   568   280], 'Color', 'w');
plot(mesh.x(:,id), bedElevation(:, id), 'k.-', 'MarkerSize', 20); hold on;
set(gca, 'XTick', [], 'YTick', []); yLims = get(gca, 'YLim');
plot(mesh.x(:, id), depth(:, id), 'k.-', xw, depth(1), '.', 'MarkerSize', 20);
plot([x1, xw, x2], [y1, eta0, eta0], 'b.-', 'Markersize', 20)
plot([x1, x1], [depth(1), y1 + 0.5], 'k--', [x2, x2], [depth(1), y1 + 0.5], 'k--',...
    [xw, xw], [depth(1), y1 + 0.5], 'k--')

ylim([depth(1)-1, yLims(2)])
t(1) = text(xw-5, depth(1) - 0.5, '${\it x}_{\omega}^{\ast}$');
t(2) = text(x1, depth(1) - 0.5, '${\it x}_{i-1/2}$');
t(3) = text(x2-7, depth(1) - 0.5, '${\it x}_{i+1/2}$');
t(4) = text(x1+3, y1-0.9, '$B_{i-1/2}$');
t(5) = text(x2-20, y2-0.3, '$B_{i+1/2}$');
t(6) = text((x1 + xw)./2, y1 - 0.5, '$\Omega_{i,a}$');
t(7) = text((x2 + xw)./2, y1 - 0.5, '$\Omega_{i,b}$');
set(t, 'Interpreter','latex', 'FontSize', 15)
set(t(6:7), 'FontSize', 20)

% refinement elements at new time level
eta1 = 0;
index = find(bedElevation(1,:).*bedElevation(end,:)<0);
id = index(1); 
x1 = mesh.x(1, id); x2 = mesh.x(end, id); 
y1 = bedElevation(1, id); y2 = bedElevation(end, id);
xw1 = x2 + (eta1 - y2)*(x1 - x2)./(y1 - y2);
yw1 = y2 + (y1 - y2)*(xw - x2)./(x1 - x2);
volume = .5*(x2 - xw)*(eta1 - y2);

eta = eta1*ones(size(mesh.x));
h = eta - bedElevation;
h(h<0) = 0;
h = Utilities.Limiter.Limiter1D.MinmodLinear(mesh,h);
eta = h + bedElevation;
eta(1, id) = bedElevation(1, id);
etaw = eta(1, id) + (eta(2, id) - eta(1, id))./(x2 - x1)*(xw - x1);
% draw
% refinement elements
figure('Position', [220   705   568   280], 'Color', 'w');
plot(mesh.x(:,id), bedElevation(:, id), 'k.-', 'MarkerSize', 20); hold on;
set(gca, 'XTick', [], 'YTick', []); yLims = get(gca, 'YLim');
plot(mesh.x(:, id), depth(:, id), 'k.-', xw, depth(1), '.', 'MarkerSize', 20);
plot([x1, xw, xw, x2], [y1, etaw, eta1, eta1], 'b.-', 'Markersize', 20)
plot([x1, x1], [depth(1), y1 + 0.5], 'k--', [x2, x2], [depth(1), y1 + 0.5], 'k--',...
    [xw, xw], [depth(1), y1 + 0.5], 'k--')

ylim([depth(1)-1, yLims(2)])
t(1) = text(xw-5, depth(1) - 0.5, '${\it x}_{\omega}^{\ast}$');
t(2) = text(x1, depth(1) - 0.5, '${\it x}_{i-1/2}$');
t(3) = text(x2-7, depth(1) - 0.5, '${\it x}_{i+1/2}$');
t(4) = text(x1+3, y1-0.9, '$B_{i-1/2}$');
t(5) = text(x2-20, y2-0.3, '$B_{i+1/2}$');
t(6) = text((x1 + xw)./2, y1 - 0.3, '$\Omega_{i,a}$');
t(7) = text((x2 + xw)./2, y1 - 0.5, '$\Omega_{i,b}$');
set(t, 'Interpreter','latex', 'FontSize', 15)
set(t(6:7), 'FontSize', 20)

% element coarsening
figure('Position', [220   705   568   280], 'Color', 'w');
plot(mesh.x(:,id), bedElevation(:, id), 'k.-', 'MarkerSize', 20); hold on;
set(gca, 'XTick', [], 'YTick', []); yLims = get(gca, 'YLim');
plot(mesh.x(:, id), depth(:, id), 'k.-', 'MarkerSize', 20);
plot([x1, x2], [y1, eta(2, id)], 'b.-', 'Markersize', 20)
plot([x1, x1], [depth(1), y1 + 0.5], 'k--', [x2, x2], [depth(1), y1 + 0.5], 'k--')

ylim([depth(1)-1, yLims(2)])
% t(1) = text(xw-5, depth(1) - 0.5, '${\it x}_{\omega}^{\ast}$');
t(1) = text(x1, depth(1) - 0.5, '${\it x}_{i-1/2}$');
t(2) = text(x2-7, depth(1) - 0.5, '${\it x}_{i+1/2}$');
t(3) = text(x1+3, y1-0.9, '$B_{i-1/2}$');
t(4) = text(x2-20, y2-0.3, '$B_{i+1/2}$');
t(6) = text((x1 + x2)./2, y1 - 0.5, '$\Omega_i$');
set(t, 'Interpreter','latex', 'FontSize', 15)
set(t(6:7), 'FontSize', 20)
set(t(6), 'FontSize', 20)