function SectionProfile_DamBreak2d
% Compare the results with exact solution on section x=0;

%% Parameters
delta = .5;
rmin  = 0+delta; 
rmax  = 1000-delta;
ne    = 1000;   % number of exact solution
np    = 60;    % number of interpolated solutions
T     = 20;

%% Get points coordinate
ye    = zeros(ne, 1);
xe    = linspace(rmin, rmax, ne)';

yp    = zeros(np, 1);
xp    = linspace(rmin, rmax, np)';

%% spicific time
timeFrac = (1/4:1/4:1);
time     = timeFrac*T;
timeStr  = cell(numel(time), 1);
for i = 1:numel(time)
    timeStr{i} = ['t= ',num2str(time(i)),' s',];
end

%% Construct postprocess class
meshtype = 'quad';
filename = {'SWE2D_DamBreakDry_quad_200.nc'};
fileID   = 1;
% create post process class for quad
PostproQuad = Utilities.PostProcess.Postprocess(filename, meshtype, 1);

meshtype = 'quad';
filename = {'SWE2D_DamBreakDry_quad_100.nc'};
PostproTri  = Utilities.PostProcess.Postprocess(filename, meshtype, 1);
lineWidth   = 2;
markerSize  = 8;
for ist = 1:numel(time)
    varname = 'h';
    extH    = DamBreakDry2d_H(xe, ye, time(ist));
    numSol  = PostproQuad.GetVarData(varname, time(ist), fileID);
    numH    = PostproQuad.Interp2D(numSol, xp, yp, fileID);
    % draw water height
    figure('Color', 'w');
    plot(xe, extH, 'k--', 'LineWidth', 1.5); hold on
    plot(xp, numH, 'ro-', 'LineWidth', lineWidth, 'MarkerSize', markerSize);
    numSol  = PostproTri.GetVarData(varname, time(ist), fileID);
    numH    = PostproTri.Interp2D(numSol, xp, yp, fileID);
    plot(xp, numH, 'bs-', 'LineWidth', lineWidth, 'MarkerSize', markerSize);
    ylabel('水位 (m)','FontSize', 16);
    xlabel('y (m)','FontSize', 16);
    title(timeStr{ist}, 'Interpreter', 'Latex', 'FontSize', 18);
    ylim([0, 11])
    t = legend('Exact', '四边形', '三角形');
    set(t, 'box', 'off', 'FontSize', 16);
    
    % draw flux
    varname = 'qx';
    extH    = DamBreakDry2d_Qx(xe, ye, time(ist));
    numSol  = PostproQuad.GetVarData(varname, time(ist), fileID);
    numH    = PostproQuad.Interp2D(numSol, xp, yp, fileID);
    figure('Color', 'w');
    plot(xe, extH, 'k--', 'LineWidth', 1.5); hold on
    plot(xp, numH, 'ro-', 'LineWidth', lineWidth, 'MarkerSize', markerSize);
    numSol  = PostproTri.GetVarData(varname, time(ist), fileID);
    numH    = PostproTri.Interp2D(numSol, xp, yp, fileID);
    plot(xp, numH, 'bs-', 'LineWidth', lineWidth, 'MarkerSize', markerSize);
    ylabel('流量 q_x (m^2/s)','FontSize', 16);
    xlabel('y (m)','FontSize', 16);
    ylim([-1, 30]);
    title(timeStr{ist}, 'Interpreter', 'Latex', 'FontSize', 18);
    t = legend('Exact', '四边形', '三角形');
    set(t, 'box', 'off', 'FontSize', 16);
    
end% for