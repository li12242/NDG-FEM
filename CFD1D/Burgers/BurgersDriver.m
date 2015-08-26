function BurgersDriver(nOrder, nElement)
% 1D Burgers¡¯ equation
% REFERENCE:
% [1] Khan A A, Lai W. Modeling Shallow Water Flows Using the 
%     Discontinuous Galerkin Method[M]. CRC Press, 2014. 50-52

% max order of base polynomial
% nOrder = 1;
% No. of element
% nElement = 100;
[Nv, VX, ~, EToV] = Utilities.MeshGen1D(-1,1,nElement);

% creat shape & mesh metric
line = StdRegions.Line(nOrder);
BC = [2, 1; 3, Nv];     % boundary condition, see MultiRegions.BCType
mesh = MultiRegions.RegionLineBC(line, EToV, VX, BC);

% Set initial conditions
initCase = 2;

% Initional condition
Q = BrugersInit(mesh, initCase);
FinalTime = 0.6;

[Q] = BurgersSolver(mesh, Q, FinalTime);
% plot(mesh.x, Q(:,:,1)); drawnow;
BurgersPostProcess(mesh, Q, initCase);
[L1, L2, Linf] = BurgersErr(mesh, Q, initCase);
fprintf('|L1 \t| L2 \t| Linf \t|\n')
fprintf('| %f \t| %f \t| %f \t|\n', L1, L2, Linf)

% [L1, L2, Linf] = BurgersErr(mesh, Q, initCase);
end% func

function [L1, L2, Linf] = BurgersErr(mesh, Q, initCase)
% analysis result
Qa = BurgersExactResult(mesh, initCase);
L1 = sum( abs(Q(:) - Qa(:) ) )./mesh.nNode;
L2 = sqrt( sum( (Q(:) - Qa(:)).^2 )/mesh.nNode );
Linf = max( abs( Q(:) - Qa(:) ) );
end% func

function Qa = BurgersExactResult(mesh, initCase)
Qa = zeros(size(mesh.x));
switch initCase
    case 1
        flag = mesh.x < 0.45;
        Qa(flag) = 1.0;
        Qa(~flag) = .5;
    case 2
        flag = mesh.x < 0.3;
        Qa(flag) = 0.5;
        flag = mesh.x > 0.6;
        Qa(flag) = 1.0;
        flag = (mesh.x >= 0.3) & (mesh.x <= 0.6);
        temp_x = mesh.x(flag);
        Qa(flag) = 5*(temp_x - 0.3)/3 + 0.5;
end% switch
end% func

function BurgersPostProcess(mesh, Q, initCase)
% Draw figure
Qa = BurgersExactResult(mesh, initCase);
plot(mesh.x(:), Qa(:),'b-'); hold on;
plot(mesh.x(:), Q(:),'ro'); 
legend('Exact', 'DGM')
end% func

function Q = BrugersInit(mesh, initCase)
flag = mesh.x < 0.0;
switch initCase
    case 1 % case 1: uL=1.0; uR=0.5;
        Q = 0.5.*ones(size(mesh.x));
        Q(flag) = 1.0;
    case 2 % case 2: uL=0.5; uR=1.0;
        Q = 1.0.*ones(size(mesh.x));
        Q(flag) = .5;
end% switch
end% func