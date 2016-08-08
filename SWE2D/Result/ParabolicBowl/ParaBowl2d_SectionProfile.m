function ParaBowl2d_SectionProfile
% Compare the results with exact solution on section x=0;

%% Parameters
g     = 9.81;
alpha = 1.6*1e-7;
w     = sqrt(8*g*alpha);
T     = 2*pi/w;
rmin  = -4000; 
rmax  =  4000;
ne    = 1000;   % number of exact solution
np    = 100;    % number of interpolated solutions


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
timeFrac = (0:1/4:1);
time     = timeFrac*T;
timeStr  = cell(numel(time), 1);
for i = 1:numel(time)
    timeStr{i} = ['t=',num2str(timeFrac(i)),'T',];
end

%% Construct postprocess class
meshtype = 'quad';
elenum   = [80, 100, 120, 160, 200];
filename = cell(numel(elenum), 1);
for i =1:numel(elenum)
    filename{i} = ['SWE2D_', meshtype, '_', num2str(elenum(i)), '.nc'];
end
fileID   = 1;
% create post process class for quad
PostproQuad = Utilities.PostProcess.Postprocess(filename, meshtype, 1);

for ist = 1:numel(time)
    varname = 'h';
    extH    = ParabolicBowlExtDepth(xe, ye, time(ist));
    numSol  = PostproQuad.GetVarData(varname, time(ist), fileID);
    numH    = PostproQuad.Interp2D(numSol, xp, yp, fileID);
    % draw water height
    figure('Color', 'w');
    plot(ye, extH+be, 'k--'); hold on
    plot(yp, numH+bp, 'r+');
    plot(ye, be, 'k');
    ylabel('$Elvation \, (m)$', 'Interpreter', 'Latex');
    xlabel('y (m)', 'Interpreter', 'Latex');
    title(timeStr{ist}, 'Interpreter', 'Latex');
    t = legend('Exact', 'RKDG');
    set(t, 'box', 'off');
    
    % draw flux
    varname = 'qx';
    extH    = ParabolicBowlExtQx(xe, ye, time(ist));
    numSol  = PostproQuad.GetVarData(varname, time(ist), fileID);
    numH    = PostproQuad.Interp2D(numSol, xp, yp, fileID);
    figure('Color', 'w');
    plot(ye, extH, 'k--'); hold on
    plot(yp, numH, 'r+');
    ylim([-1.25, 1.25]);
    ylabel('$\rm{Discharge} \, q_x \, (m^2/s)$', 'Interpreter', 'Latex');
    xlabel('y (m)', 'Interpreter', 'Latex');
    title(timeStr{ist}, 'Interpreter', 'Latex');
    t = legend('Exact', 'RKDG');
    set(t, 'box', 'off');
    
    % draw flux
    varname = 'qy';
    extH    = ParabolicBowlExtQx(xe, ye, time(ist));
    numSol  = PostproQuad.GetVarData(varname, time(ist), fileID);
    numH    = PostproQuad.Interp2D(numSol, xp, yp, fileID);
    figure('Color', 'w');
    plot(ye, extH, 'k--'); hold on
    plot(yp, numH, 'r+');
    ylim([-1.25, 1.25]);
    ylabel('$\rm{Discharge} \, q_x \, (m^2/s)$', 'Interpreter', 'Latex');
    xlabel('y (m)', 'Interpreter', 'Latex');
    title(timeStr{ist}, 'Interpreter', 'Latex');
    t = legend('Exact', 'RKDG');
    set(t, 'box', 'off');
end% for