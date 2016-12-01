function SectionProfile_ParaBowl2d
% Compare the results with exact solution on section x=0;

%% Parameters
g     = 9.81;
alpha = 1.6*1e-7;
w     = sqrt(8*g*alpha);
T     = 2*pi/w;
delta = 0.5;
rmin  = -4000+delta; 
rmax  =  4000-delta;
ne    = 100;   % number of exact solution
np    = 50;    % number of interpolated solutions

%% Get points coordinate
xe    = zeros(ne, 1);
ye    = linspace(rmin, rmax, ne)';
re    = (xe.^2 + ye.^2);
be    = alpha*re;

xp    = zeros(np, 1);  
yp    = linspace(rmin, rmax, np)';
rp    = (xp.^2 + yp.^2);
bp    = alpha*rp;

%% spicific time
timeFrac = linspace(eps, 1-eps, 4);
% timeFrac = (1/4:1/4:3/4);
time     = timeFrac*T;
timeStr  = cell(numel(time), 1);
for i = 1:numel(time)
    timeStr{i} = ['t=',num2str(timeFrac(i)),'T',];
end

%% Construct postprocess class
meshtype = 'quad';
filename = {'SWE2D_ParabolicBowl_quad_60.nc'};
Postpro1 = Utilities.PostProcess.Postprocess(filename, meshtype, 1);

meshtype = 'tri';
filename = {'SWE2D_ParabolicBowl_tri_60.nc'};
Postpro2 = Utilities.PostProcess.Postprocess(filename, meshtype, 1);

lineWidth = 2;
markerSize = 7;
for ist = 1:numel(time)
    varname = 'h';
    extH    = ParabolicBowl2d_ExtDepth(xe, ye, time(ist));
    numH1 = InterpSolution(Postpro1, varname, time(ist), xp, yp);
    numH2 = InterpSolution(Postpro2, varname, time(ist), xp, yp);
    % draw water height
    figure('Color', 'w');
    plot(ye, extH+be, 'k--', 'LineWidth',lineWidth); hold on
    plot(yp, numH1+bp, 'ro-', 'MarkerSize', markerSize,...
        'LineWidth',lineWidth, 'MarkerFaceColor', 'r');
    plot(yp, numH2+bp, 'b^-', 'MarkerSize', markerSize,...
        'LineWidth',lineWidth, 'MarkerFaceColor', 'b');
    plot(ye, be, 'k', 'LineWidth',lineWidth);
    ylabel('\eta (m)', 'FontSize', 14);
    xlabel('y (m)', 'FontSize', 14);
    title(timeStr{ist}, 'FontSize', 14);
    legend({'Exact', 'Quad', 'Tri'}, 'box', 'off','FontSize', 14);
    
    % draw flux qx
    varname = 'qx';
    extH    = ParabolicBowl2d_ExtQx(xe, ye, time(ist));
    numH1 = InterpSolution(Postpro1, varname, time(ist), xp, yp);
    numH2 = InterpSolution(Postpro2, varname, time(ist), xp, yp);
    figure('Color', 'w');
    plot(ye, extH, 'k--', 'LineWidth',lineWidth); hold on
    plot(yp, numH1, 'ro-', 'MarkerSize', markerSize, 'LineWidth',lineWidth,...
        'MarkerFaceColor', 'r');
    plot(yp, numH2, 'b^-', 'MarkerSize', markerSize, 'LineWidth',lineWidth,...
        'MarkerFaceColor', 'b');
    ylim([-1.25, 1.25]);
    ylabel('q_x (m^2/s)', 'FontSize', 14);
    xlabel('y (m)', 'FontSize', 14);
    title(timeStr{ist}, 'FontSize', 14);
    legend({'Exact', 'Quad', 'Tri'}, 'box', 'off','FontSize', 14);
    
    % draw flux qy
    varname = 'qy';
    extH    = ParabolicBowl2d_ExtQy(xe, ye, time(ist));
    numH1 = InterpSolution(Postpro1, varname, time(ist), xp, yp);
    numH2 = InterpSolution(Postpro2, varname, time(ist), xp, yp);
    figure('Color', 'w');
    plot(ye, extH, 'k--', 'LineWidth',lineWidth); hold on
    plot(yp, numH1, 'ro-', 'MarkerSize', markerSize, 'LineWidth',lineWidth,...
        'MarkerFaceColor', 'r');
    plot(yp, numH2, 'bs-', 'MarkerSize', markerSize, 'LineWidth',lineWidth,...
        'MarkerFaceColor', 'b');
    ylim([-1.25, 1.25]);
    ylabel('q_y (m^2/s)', 'FontSize', 14);
    xlabel('y (m)', 'FontSize', 14);
    title(timeStr{ist}, 'FontSize', 14);
    legend({'Exact', 'Quad', 'Tri'}, 'box', 'off', 'FontSize', 14);
end% for
end% func

function interpSol = InterpSolution(Postpro, varname, time, xp, yp)
fileID  = 1;
numSol  = Postpro.GetVarData(varname, time, fileID);
interpSol= Postpro.Interp2D(numSol, xp, yp, fileID);
end