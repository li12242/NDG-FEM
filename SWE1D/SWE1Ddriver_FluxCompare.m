% Driver script for solving the 1D shallow water equation

% Order of polymomials used for approximation 
N = 3; 

nElement = 60;
assert(mod(nElement,2) == 0)
% Generate simple mesh
[Nv, VX, K, EToV] = Utilities.MeshGen1D(0.0,1000,nElement);

% Initialize solver and construct grid and metric
point = StdRegions.Point();
line = StdRegions.Line(N, point);
mesh = MultiRegions.RegionLine(line, EToV, VX);

% Solve Problem
FinalTime = 30;

%
figure('Color', 'w');
hold on;
lineColor = {'r','b', 'g', 'c'};

fluxFunc{1} = @SWE1DRHS1D;
fluxFunc{2} = @SWE1DRHS_LaxFriedrich;
fluxFunc{3} = @SWE1DRHS_RUS;
fluxFunc{4} = @SWE1DRHS_HLL;
% tic
for i = 1:4
    % Set initial conditions
    Q = zeros([size(mesh.x),2]);
    Q(:,:,1) = 2.1.*ones(size(mesh.x));
    Q(:,1:nElement/2, 1) = 10;
    Q(:,:,2) = zeros(size(mesh.x));
    
    [Q] = SWE1D(mesh,Q,FinalTime, fluxFunc{i});
    plot(mesh.x, Q(:,:,1), [lineColor{i}, 'o']); drawnow;
end
% toc

pix_left = 126; pix_right = 706;
pix_top = 72; pix_bottom = 601;
ratio.x = 1000/(pix_right - pix_left); ratio.y = 10./(pix_bottom-pix_top);
P(1,:) = [126-pix_left, pix_bottom-72].*[ratio.x, ratio.y];
P(2,:) = [240-pix_left, pix_bottom-72].*[ratio.x, ratio.y];
P(3,:) = [387-pix_left, pix_bottom-332].*[ratio.x, ratio.y];
P(4,:) = [578-pix_left, pix_bottom-332].*[ratio.x, ratio.y];
P(5,:) = [578-pix_left, pix_bottom-493].*[ratio.x, ratio.y];
P(6,:) = [706-pix_left, pix_bottom-493].*[ratio.x, ratio.y];
plot(P(:,1), P(:,2), 'ko-');
set(gca, 'YLim', [0,11]);
% legend('Global Lax-Friedrich', 'Local L-F', 'Rusanov','Theoretial')