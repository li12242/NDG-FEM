function var = Convection2DSetUp
% 2D convection problem
% dc/dt + udc/dx + vdu/dy = 0

N = 2;
% mesh = quadSolver(N);
mesh = triSolver(N);
var = ConvectionInit(mesh);

Speed = [1,0]; % speed of domain, [u, v]
FinalTime = 0.25;

var = Convection2DSolver(mesh, var, FinalTime, Speed);
postprocess(mesh, var)
end% func

function var = ConvectionInit(mesh)
% var = ones(size(mesh.x));
% var = mesh.x;
% xc = mean(mesh.x);
% left = xc < 0.5; right = xc > 0.5;
var = sin(pi*mesh.x);%.*sin(2*pi*mesh.y);
% var = zeros(size(mesh.x));
% var(:,left) = 1; var(:,right) = 0;
end% func

function mesh = triSolver(N)

% read triangle mesh
% [EToV, VX, VY, EToR, BC] = Utilities.Mesh.MeshReaderTriangle('Convection2D/mesh/triangle');
M = 20;
np = M+1;
[VX,VY,EToV] = Lmesh(np, -1, 1, 1);

tri = StdRegions.Triangle(N);
mesh = MultiRegions.RegionTri(tri, EToV, VX, VY);
end% func

function mesh = quadSolver(N)

% read triangle mesh
% [EToV, VX, VY, EToR, BC] = Utilities.Mesh.MeshReaderQuad('Convection2D/mesh/quad');

% uniform mesh
[EToV, VX, VY] = uniformRetangle;

temp = EToV(:, 3); 
EToV(:, 3) = EToV(:, 4);
EToV(:, 4) = temp;

quad = StdRegions.Quad(N);
mesh = MultiRegions.RegionQuad(quad, EToV, VX, VY);

end% func

function [X,Y,E] = Lmesh(n,start_coor,end_coor,flag)
% 绘制规则三角形地形
%
% DESCRIPTION
%   区域范围 [start_coor, end_coor]x[start_coor, end_coor];
%   每边单元个数为n，节点个数为n+1
%
%   生成网格节点编号顺序为:
%       由(0, 1)起始，沿x轴方向循环
%       1   2   3   4   5   ....            n  n+1
%       -----------------------------------------
%       |   |   |   |   |   |   |   |   |   |   |
%
% INPUT
%    n      = 每边单元个数
%    start_coor = 起始坐标 
%    end_coor   = 结束坐标
%    flag   = [flag=0, 斜边方向为"\"；flag=1，斜边方向为"/"]
%
% OUTPUT
%    X = x coordinate
%    Y = y coordinate
%    E = points of element (逆时针)
%
% EXAMPLE USAGE
%    [X,Y,E] = Lmesh(21, 1)
%
% Author(s)
%    li12242 Tianjin University
%
% Reversion
%    v1.0 2014-12-8
%========================================================================== 

% info of point
x = linspace(start_coor, end_coor, n+1); y = linspace(end_coor, start_coor, n+1);
x = repmat(x, n+1, 1); y = repmat(y, n+1, 1);
x = x';
X = x(:); Y = y(:);
% 单元组成
E = zeros(2*n^2,3);
for i = 1:n %每层单元循环
    upLayerNum = (n+1)*(i-1)+1:(n+1)*(i-1)+n+1;
    downLayerNum = upLayerNum + (n+1);

    if flag     % 上层单元在先
        E(2*n*(i-1)+1:2*n*(i-1)+n,:) = [downLayerNum(1:n)', upLayerNum(2:n+1)', upLayerNum(1:n)'];
        E(2*n*(i-1)+n+1:2*n*i,:) = [downLayerNum(1:n)',downLayerNum(2:n+1)',upLayerNum(2:n+1)'];
        flag=~flag;
    else
        E(2*n*(i-1)+1:2*n*(i-1)+n,:) = [downLayerNum(2:n+1)', upLayerNum(2:n+1)', upLayerNum(1:n)'];
        E(2*n*(i-1)+n+1:2*n*i,:) = [downLayerNum(1:n)',downLayerNum(2:n+1)',upLayerNum(1:n)'];
        flag=~flag;
    end

end
end

function [EToV, VX, VY] = uniformRetangle(np)
m = np - 1; % elements on each row

VX = 1:np; VX = repmat(VX, 1, np);
VX = linspace(-1, 1, np);
dx = [-1, 0, -1, 0, -1];
for i = 1:np
    VX((i-1)*np+1: i*np) = VX((i-1)*np+1: i*np) + dx(i);
end
VY = linspace(1, -1, np)'; VY = repmat(VY, 1, np)'; VY = VY(:);

EToV = zeros(m^2, 4);
for irow = 1:m
    for icol = 1:m
        index = (irow-1)*m + icol;
        EToV(index, :) = [np*(irow)+icol, np*(irow)+icol+1, ...
            np*(irow-1)+icol, np*(irow-1)+icol+1];
%         EToV(index, :) = [np*(irow-1)+icol, np*(irow-1)+icol+1,...
%             np*(irow)+icol, np*(irow)+icol+1];
    end
end

end% func

function postprocess(mesh, var)
plot3(mesh.x(mesh.vmapP),mesh.y(mesh.vmapP), var(mesh.vmapM))
end% func

