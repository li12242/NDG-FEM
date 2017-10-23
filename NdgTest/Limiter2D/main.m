function main
close all
% test of 2D slope limiter
test = 1; 

switch test
    case 1 % quad
        M = 2; N = 1;
        mesh = quadSolver(N, M);
    case 2 % triangle 
        M = 2; N = 1;
        mesh = triSolver(N, M);
end

var = ConvectionInit(mesh);

% variable distribution
figure('Position', [680,1,560,975]); subplot(2,1,1); hold on;
plot3(mesh.x(mesh.vmapP),mesh.y(mesh.vmapP), var(mesh.vmapM), 'k-');
zlim([-1, 1]);
view([20, 50]);

% discontinuity detector
% var = Utilities.Limiter.Limiter2D.HWENO2d(mesh, var);
ind    = Utilities.Limiter.Limiter2D.TVB_detector2d(mesh, var, 0);
% varlim = Utilities.Limiter.Limiter2D.VB2d(mesh, var);
varlim = Utilities.Limiter.Limiter2D.TVB_tri2d(mesh, var, 0.2);

subplot(2,1,2); hold on;
plot3(mesh.x(mesh.vmapP),mesh.y(mesh.vmapP), varlim(mesh.vmapM), 'k-');
zlim([-1.5, 1.5]);
view([20, 50]);
end% func

function var = ConvectionInit(mesh)
shape = mesh.Shape;
hmean = 0;
w     = sum(shape.M);
area  = (w*mesh.J);
xmean = (w*(mesh.J.*mesh.x))./area;
ymean = (w*(mesh.J.*mesh.y))./area;

xmean = repmat(xmean, shape.nNode, 1);
ymean = repmat(ymean, shape.nNode, 1);

gra   = [1, 0];
var = hmean+(mesh.x - xmean).*gra(1) + (mesh.y - ymean).*gra(2);
end

function var = Square(mesh, x0, y0, r0)
var = zeros(size(mesh.x));
% slotted cylinder
ind = (abs(mesh.x-x0)<r0 & abs(mesh.y-y0)<r0);
var(ind) = 1.0;
end% func

function var = Cylinder(mesh, x0, y0, r0)
var = zeros(size(mesh.x));
% slotted cylinder
r2  = sqrt((mesh.x-x0).^2+(mesh.y-y0).^2)./r0;
ind = ( r2<=1.0);
var(ind) = 1.0;
end% func

function g = G(x, y, beta, x0, y0)
r2 = (x - x0).^2 + (y - y0).^2;
g = exp( -beta*r2 );
end

function f = F(x, y, alpha, a, b)
r2 = (x - a).^2 + (y - b).^2;
f = sqrt( max(1 - alpha*r2, 0) );
end% func

function mesh = triSolver(N, M)

[EToV, VX,VY] = Utilities.Mesh.MeshGenTriangle2D(M+1, M+1, -1, 1, -1, 1, 0);
tri = StdRegions.Triangle(N);
mesh = MultiRegions.RegionTri(tri, EToV, VX, VY);
end% func

function mesh = quadSolver(N, M)
% uniform mesh
[EToV, VX, VY] = Utilities.Mesh.MeshGenRectangle2D(M+1, M+1, -1, 1, -1, 1);

quad = StdRegions.Quad(N);
mesh = MultiRegions.RegionQuad(quad, EToV, VX, VY);

end% func