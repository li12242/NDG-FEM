function [h, q] = SWE_PositivePreserving1d(phys, h, q)
% Slope limiter and Positivity-preserving operator
mesh = phys.mesh;
%% slope limiter on water level and discharge

fp = fopen('test.txt', 'a+');
fmtstr = '%f %f\n';
fprintf(fp, '\nh_ = \n');
fprintf(fp, fmtstr, h);
fprintf(fp, '\nq_ = \n');
fprintf(fp, fmtstr, q);

% the slope limiter act on the wet cells
q = Utilities.Limiter.Limiter1D.MinmodLinear(mesh, q);
h = Utilities.Limiter.Limiter1D.MinmodLinear(mesh, h); 

fprintf(fp, '\nh_1 = \n');
fprintf(fp, fmtstr, h);
fprintf(fp, '\nq_1 = \n');
fprintf(fp, fmtstr, q);

%% positive preserving operator
%shape  = mesh.Shape;
hc = mean(h); qc = mean(q);
% [h, q] = SWE_Mex_PositivePreserving1d(h, q, shape.M, mesh.J, phys.minDepth);
[h, q] = ppreserve(phys.minDepth, h, q, hc, qc);
% [h, q] = PositiveOperator(mesh, h, q);

fprintf(fp, '\nh_2 = \n');
fprintf(fp, fmtstr, h);
fprintf(fp, '\nq_2 = \n');
fprintf(fp, fmtstr, q);
fprintf(fp, '\n======================\n');
fclose(fp);
end% func

function [h, q] = PositiveOperator(mesh, h, q)
% positive operator
% reference from Xing (2010); Zhang (2010)

ksi    = 0.0;
hmean  = SWE_CellMean1d(mesh, h);
Np     = mesh.Shape.nNode;
% correct mean water less than hDelta
dis    = (hmean <= ksi);
h(:, dis) = h(:, dis) + ones(Np, 1)*(ksi - hmean(dis));
q(:, dis) = 0;

hmean = SWE_CellMean1d(mesh, h);
qmean = SWE_CellMean1d(mesh, q);
% positive operator
hmin = min(h);
theta = min( (hmean - ksi)./(hmean - hmin), 1);
h = (ones(Np, 1)*theta).*(h - ones(Np, 1)*hmean) + ones(Np, 1)*hmean;
q = (ones(Np, 1)*theta).*(q - ones(Np, 1)*qmean) + ones(Np, 1)*qmean;
end% func