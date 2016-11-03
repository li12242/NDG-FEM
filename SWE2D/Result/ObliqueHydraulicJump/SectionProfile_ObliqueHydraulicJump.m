function SectionProfile_ObliqueHydraulicJump
%SECTIONPROFILE_OBLIQUEHYDRAULICJUMP Summary of this function goes here
%   Detailed explanation goes here

np = 100;
y0 = 20;
EPS = 1e-10;
xp = linspace(eps, 40-EPS, np);
yp = ones(1, np)*y0;
hp = ExactH(xp, yp);

figure('Color', 'w');
plot(xp, hp, 'k-', 'LineWidth', 2); hold on;

% Draw result
meshtype = 'tri';
filename = ['SWE2D_ObliqueHydraulicJump_VBVA_',meshtype,'.nc'];
p0 = DrawProfile(y0, filename, meshtype, 'r--');

filename = ['SWE2D_ObliqueHydraulicJump_VBJK_',meshtype,'.nc'];
% filename = 'SWE2D_ObliqueHydraulicJump_quad_300.nc';
p1 = DrawProfile(y0, filename, meshtype, 'b-.');

filename = ['SWE2D_ObliqueHydraulicJump_VBHWENO_',meshtype,'.nc'];
p2 = DrawProfile(y0, filename, meshtype, 'g:');

legend({'Exact','VA','JK','QS'}, ...
    'Location', 'NorthWest', 'box', 'off', 'FontSize', 16)

ylabel('$h (m)$','Interpreter', 'Latex', 'FontSize', 16);
xlabel('$x (m)$','Interpreter', 'Latex', 'FontSize', 16);
box on; grid on;
title(['$y=',num2str(y0),' m$'],'Interpreter', 'Latex', 'FontSize', 18);
end

function hp = ExactH(xp, yp)
k = -tan(30/180*pi);
hp = ones(size(xp));
flag = ((xp - 10)*k + 30)<yp;
hp(flag) = 1.5049;
end

function p = DrawProfile(y0, filename, meshType, lineType)
% Parameter
lineWidth   = 2.5;
markerSize  = 8;
finalTime   = 14;
fileID      = 1;
np = 50;
xp = linspace(eps, 40-eps, np);
yp = ones(1, np)*y0;
varname = 'h';

strcname = {filename};
Postpro  = Utilities.PostProcess.Postprocess(strcname, meshType, fileID);
numSol  = Postpro.GetVarData(varname, finalTime, fileID);
numH    = Postpro.Interp2D(numSol, xp, yp, fileID);
p = plot(xp, numH, lineType, 'LineWidth', lineWidth, ...
    'MarkerSize', markerSize);
end% func

