function drawMaxRunup( obj )

bctype = [...
    NdgEdgeType.ZeroGrad, ...
    NdgEdgeType.ZeroGrad, ...
    NdgEdgeType.ClampedDepth, ...
    NdgEdgeType.ZeroGrad];

M = 150;
N = 1;
mesh = makeUniformQuadMesh(N, [0, 25], [0, 30], M, M, bctype);
%mesh = obj.meshUnion;

casename = 'ConicalLand2d_150_1';
conpos = NdgPostProcess( mesh, casename );

dryFlag = ones( mesh.cell.Np, mesh.K );
threshold = 2.3e-3;
for t = 1:conpos.Nt
    [ fphys ] = conpos.accessOutputResultAtStepNum( t );
    dep = fphys{1}(:,:,1);
    dryFlag( dep > threshold ) = 0;
end

maxdata = [10.9361,13.8192
11.0819,13.0178
11.0929,14.5888
11.4955,15.2634
11.5175,12.3486
12.1546,11.8621
12.1694,15.7422
12.9552,11.6846
12.9633,15.9117
13.7744,11.8162
13.7765,15.7973
14.4765,15.3423
14.4924,12.2823
14.9255,13.8907
14.9319,13.8087
14.932,13.733
14.9506,13.9916
14.9563,14.6352
14.9636,13.632
14.9643,12.9822
15.052,13.4489
15.0576,14.1934
15.0649,13.2407
15.07,14.3701];

dx = - 0.42; 
dy = 1.2;

dryFlag( dryFlag < 1 ) = nan;
plot( maxdata(:, 1)+dx, maxdata(:, 2)+dy, 'rx', 'MarkerSize', 10); hold on;

plot3( mesh.x(:), mesh.y(:), dryFlag(:), 'b.', 'MarkerSize', 10)
box on;
view([0, 90]);
axis equal;

legend({'Measured max run-up point', 'Computed dry mesh nodes with $p=1$'}, ...
    'Interpreter', 'latex', 'FontSize', 16, 'box', 'off', ...
    'Location', 'northoutside');
xlabel('$x$ (m)', 'Interpreter', 'latex','FontSize', 16);
ylabel('$y$ (m)', 'Interpreter', 'latex','FontSize', 16);
end

