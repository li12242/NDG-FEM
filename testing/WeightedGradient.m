quad = StdRegions.Quad(2);
VX = [0, 1, 1, 0]';
VY = [0, 0, 1, 1]';
EToV = [1,2,3,4];
mesh = MultiRegions.RegionQuad(quad, EToV, VX, VY);
h = zeros(size(mesh.x));
h = mesh.x + rand(size(mesh.x));

gra = @(x, y, h) ...
    ([x(2)-x(1), y(2)-y(1); x(3)-x(1), y(3)-y(1)])\[h(2) - h(1); h(3)-h(1)];

vertlist = quad.getVertexNodeList;
w  = sum(quad.M)';
area = sum(w.*mesh.J);
xc = sum(w.*mesh.J.*mesh.x)/area;
yc = sum(w.*mesh.J.*mesh.y)/area;
hc = sum(w.*mesh.J.*h)/area;

dh = zeros(2,4);
for i = 1:4;
    l1 = i; l2 = mod(i, 4)+1;
    x1 = mesh.x(vertlist(l1)); x2 = mesh.x(vertlist(l2));
    y1 = mesh.y(vertlist(l1)); y2 = mesh.y(vertlist(l2));
    h1 = h(vertlist(l1));      h2 = h(vertlist(l2));
    dh(:, i) = gra([x1, x2, xc], [y1, y2, yc], [h1, h2, hc]);
end% for

g = sum(dh.^2);
epse = 1e-10;
% frac = g(2)*g(3)*g(4)+g(1)*g(3)*g(4)+g(2)*g(1)*g(4)+g(2)*g(3)*g(1);
frac = sum(g.^3);
wg = ([g(2)*g(3)*g(4), g(1)*g(3)*g(4), g(2)*g(1)*g(4), ...
    g(2)*g(3)*g(1)] + epse)/(frac+3*epse);

dhlim = sum(repmat(wg, 2, 1).*dh, 2);

hlim = (mesh.x-xc).*dhlim(1) + (mesh.y-yc).*dhlim(2) + hc;

figure('Color', 'w')
p1 = plot3(mesh.x(mesh.vmapM), mesh.y(mesh.vmapM), h(mesh.vmapM)); hold on;
p2 = plot3(mesh.x(mesh.vmapM), mesh.y(mesh.vmapM), hlim(mesh.vmapM));
legend([p1, p2], {'original', 'limited'}, 'Location', 'NorthWest', 'box', 'off')