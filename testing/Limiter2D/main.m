function main
% test of 2D slope limiter
test = 1; 

switch test
    case 1 % quad
        M = 40; N = 2;
        mesh = quadSolver(N, M);
    case 2 % triangle 
        M = 1; N = 2;
        mesh = triSolver(N, M);
end

var = ConvectionInit(mesh);
u = zeros(size(mesh.x)); v = ones(size(mesh.x));

% variable distribution
figure;
plot3(mesh.x(mesh.vmapP),mesh.y(mesh.vmapP), var(mesh.vmapM));

% discontinuity detector
[~, I] = Utilities.Limiter.Limiter2D.DisDetector(mesh, var, u, v);

figure
xc = mean(mesh.x); yc = mean(mesh.y);
plot3(xc, yc, I, 'r.');
end% func


function var = ConvectionInit(mesh)
% initial condition
% refer to Krivodonova (2007) for more details
var = zeros(size(mesh.x));

% up left
delta = 0.15;
beta = log(2)./delta^2;

temp = G(mesh.x, mesh.y, beta, -.5, 0.5);

x1 = -0.8; x2 = -0.2; y1 = 0.2; y2 = 0.8;
xm = (x1 + x2)./2; ym = (y1 + y2)./2;
R = (x1 - xm).^2 + (y1 - ym).^2;
ind = ( (mesh.x - xm).^2 + (mesh.y - ym).^2 ) < R;
var(ind) = temp(ind);

% up right
x1 = 0.25; x2 = 0.75; y1 = 0.25; y2 = 0.75;
xc = mean(mesh.x); yc = mean(mesh.y);
ind = (xc > x1) & (xc < x2) ...
    & ( yc > y1 ) & (yc < y2);
var(:, ind) = 1;

% down left
x1 = -0.25; x2 = -0.75; y1 = -0.75; y2 = -0.25;
xm = (x1 + x2)./2; ym = (y1 + y2)./2;
R = (x1 - xm).^2 + (y1 - ym).^2;
ind = ( (mesh.x - xm).^2 + (mesh.y - ym).^2 ) < R;

temp = 1 - sqrt( ((mesh.x - xm).^2 + (mesh.y - ym).^2)./R );
var(ind) = temp(ind);

% down right
x1 = 0.25; x2 = 0.75; y1 = -0.75; y2 = -0.25;
xm = (x1 + x2)./2; ym = (y1 + y2)./2;
R = (x1 - xm).^2 + (y1 - ym).^2;
ind = ( (mesh.x - xm).^2 + (mesh.y - ym).^2 ) < R;

alpha = 10;
xm = 0.5; ym = -0.5;
temp = F(mesh.x, mesh.y, alpha, xm, ym);
var(ind) = temp(ind);
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

[VX,VY,EToV] = Utilities.Mesh.MeshGenTriangle2D(M, -1, 1, 1);
tri = StdRegions.Triangle(N);
mesh = MultiRegions.RegionTri(tri, EToV, VX, VY);
end% func

function mesh = quadSolver(N, M)
% uniform mesh
[EToV, VX, VY] = Utilities.Mesh.MeshGenRectangle2D(M+1, -1, 1);

quad = StdRegions.Quad(N);
mesh = MultiRegions.RegionQuad(quad, EToV, VX, VY);

end% func



