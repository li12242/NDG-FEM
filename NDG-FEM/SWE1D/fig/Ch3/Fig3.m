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

% partially wet transition element
eta0 = -2;
index = find(bedElevation(1,:).*bedElevation(end,:)<0);
id = index(1); 
x1 = mesh.x(1, id); x2 = mesh.x(end, id); 
y1 = bedElevation(1, id); y2 = bedElevation(end, id);
xw = x2 + (eta0 - y2)*(x1 - x2)./(y1 - y2);
yw = y2 + (y1 - y2)*(xw - x2)./(x1 - x2);
volume = .5*(x2 - xw)*(eta0 - y2);
hmean = volume./(x2 - x1);
etaw = 0.5*(y1 + y2) + hmean;

eta = eta0*ones(size(mesh.x));
depth = (eta0 - 3)*ones(size(mesh.x));
% draw figure
close all; figure('Position', [220   705   568   280], 'Color', 'w');
plot([x1, xw, x2], [y1, yw, eta0], 'b')
plot(mesh.x(:,id), bedElevation(:, id), 'k.-', 'MarkerSize', 20)
hold on; set(gca, 'XTick', [], 'YTick', []); yLims = get(gca, 'YLim');
plot(mesh.x(:, id), depth(:, id), 'k.-', xw, depth(1), '.', 'MarkerSize', 20);

xlist = [xw, x2, x2]'; ylist = [yw, y2, eta0]';
fill(xlist, ylist, 'b')

plot([x1, (x1 + x2)./2, x2], [y1, etaw, y2+2*hmean], 'r.-', 'MarkerSize', 20);

ylim([depth(1)-1, yLims(2)])

t(1) = text(xw-5, depth(1) - 0.5, '${\it x}_{\omega}^{\ast}$');
t(2) = text(x1, depth(1) - 0.5, '${\it x}_{i-1/2}$');
t(3) = text(x2-7, depth(1) - 0.5, '${\it x}_{i+1/2}$');
t(4) = text(xw + 40, eta0+0.5, '$\Delta x \cdot \bar{h}_i$');
set(t(4), 'BackgroundColor', 'w')
t(5) = text(x1+5, y1, '$B_{i-1/2}$');
t(6) = text(x2-10, y2-0.3, '$B_{i+1/2}$');
t(7) = text((x1 + x2)./2 - 10, etaw+0.5, '$B_i + \bar{h}_i$', 'Interpreter','latex');
t(8) = text(x2 + 10, y2 + 0.4 , '$2\bar{h}_i$' );
t(9) = text(x2 + 7, (y2 + 1.2), '$h_{i+1/2}^+$' );
set(t, 'Interpreter','latex', 'FontSize', 15)
xb = get(gca, 'XLim');
xlim([xb(1), xb(2)+20])
plot([x2+4, x2+8], [eta0, eta0], 'k', [x2+4, x2+8], [y2, y2], 'k')
plot([x2+11, x2+15], [y2, y2], 'k', [x2+11, x2+15], [y2+2*hmean, y2+2*hmean], 'k')


% fully flooded transition element
eta1 = 2;

eta = eta1*ones(size(mesh.x));
h1 = eta(1, id) - bedElevation(1, id); h2 = eta(2, id) - bedElevation(2, id);
hmean = 0.5*(h1 + h2);
etaw = hmean + mean(bedElevation(:, id));
depth = (eta0 - 3)*ones(size(mesh.x));
y1 = bedElevation(1, id); y2 = bedElevation(end, id);

% draw figure
figure('Position', [220   705   568   280], 'Color', 'w');
plot(mesh.x(:, id), eta(:, id), 'b', mesh.x(:,id), bedElevation(:, id), 'k')
hold on; set(gca, 'XTick', [], 'YTick', []);
plot(mesh.x(:, id), depth(:, id), 'k.-', 'MarkerSize', 20);

xlist = [x1, x1, x2, x2]'; ylist = [eta1, y1, y2, eta1]';
fill(xlist, ylist, 'b')
plot([x1, (x1 + x2)./2, x2], [y1, (y1 + y2)./2+hmean, y2+2*hmean], 'r.-', 'MarkerSize', 20);

ylim([depth(1)-1, eta1+1])

t(1) = text(x1, depth(1) - 0.5, '${\it x}_{i-1/2}$');
t(2) = text(x2-7, depth(1) - 0.5, '${\it x}_{i+1/2}$');
t(3) = text((x1+x2)./2 - 10, eta1 + 0.5, '$B_i + \bar{h}_i$');
t(4) = text((x1+x2)./2 - 5, eta1 - 1.7, '$\Delta x \cdot \bar{h}_i$');
set(t(4), 'BackgroundColor', 'w')
t(5) = text(x1, y1-0.7, '$B_{i-1/2}$');
t(6) = text(x2-10, y2-0.3, '$B_{i+1/2}$');
t(8) = text(x2 + 10, y2 + 0.4 , '$2\bar{h}_i$' );
t(9) = text(x2 + 7, (y2 + 1.2), '$h_{i+1/2}^+$' );
set(t, 'Interpreter','latex', 'FontSize', 15)

xb = get(gca, 'XLim');
xlim([xb(1), xb(2)+30])

dx1 = 4; dx2 = 25.5;
plot([x2+dx1, x2+dx1+4], [eta1, eta1], 'k', [x2+dx1, x2+dx1+4], [y2, y2], 'k')
plot([x2+dx2, x2+dx2+4], [y2, y2], 'k', [x2+dx2, x2+dx2+4], [y2+2*hmean, y2+2*hmean], 'k')
