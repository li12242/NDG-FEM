function SectionProfile_ObliqueHydraulicJump
%SECTIONPROFILE_OBLIQUEHYDRAULICJUMP Summary of this function goes here
%   Detailed explanation goes here

np = 50;
y0 = 20;
xp = linspace(eps, 40-eps, np);
yp = ones(1, np)*y0;
hp = ExactH(xp, yp);

figure('Color', 'w');
plot(xp, hp, 'k--', 'LineWidth', 1.5); hold on;

% Draw result
filename = 'SWE2D_ObliqueHydraulicJump_tri_80.nc';
p0 = DrawProfile(y0, filename, 'tri', 'r^-');

filename = 'ObliqueHydraulicJump_VBHWENO_quad.nc';
p2 = DrawProfile(y0, filename, 'quad', 'go-');

filename = 'ObliqueHydraulicJump_VBHWENO_tri.nc';
p1 = DrawProfile(y0, filename, 'tri', 'bs-');

legend({'Exact','四边形','三角形(稀疏)','三角形(加密)'}, ...
    'Location', 'NorthWest', 'box', 'off', 'FontSize', 16)

ylabel('水位 (m)','FontSize', 16);
xlabel('y (m)','FontSize', 16);
box on; grid on;
title(['y=',num2str(y0),' m'],'Interpreter', 'Latex', 'FontSize', 18);
end

function hp = ExactH(xp, yp)
k = -tan(30/180*pi);
hp = ones(size(xp));
flag = ((xp - 10)*k + 30)<yp;
hp(flag) = 1.5049;
end

function p = DrawProfile(y0, filename, meshType, lineType)
% Parameter
lineWidth   = 2;
markerSize  = 8;
finalTime   = 15;
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

