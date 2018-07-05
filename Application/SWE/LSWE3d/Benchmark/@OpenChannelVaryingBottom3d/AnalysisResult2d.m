function AnalysisResult2d( obj )
%ANALYSISRESULT Summary of this function goes here
%   Detailed explanation goes here

OutstepNum = obj.outputFile.outputStep;

% analysis 2D water level results
cellId2 = [ 1, 6, 11, 20 ];
nodeId2 = [ 1, 1, 1, 2 ];
ind2 = sub2ind( [obj.mesh2d.cell.Np, obj.mesh2d.K], nodeId2, cellId2 );
Nx2 = numel(ind2);
x2 = obj.mesh2d.x(ind2);

etaNum = zeros( Nx2, OutstepNum );
etaExt = ExactWaterLevel( obj, x2, obj.outputFile.outputTime );

% analysis 3d horizontal velocity
% Nz = obj.mesh3d.Nz;
% cellId3 = (cellId2 - 1) * Nz + 1;
% nodeId3 = [ 1, 1, 1, 2 ];
% ind3 = sub2ind( [obj.mesh3d.cell.Np, obj.mesh3d.K], nodeId3, cellId3 );

% uNum = zeros( Nx2, OutstepNum );
% uExt = zeros( Nx2, OutstepNum );


for n = 1 : OutstepNum
    [ fphys2d, ~ ] = obj.outputFile.readOutputResult(n);  
    etaNum( :, n ) = fphys2d( ind2 );
    
%     uNum( :, n ) = fphys3d( ind3 );
%     uExt( :, n ) = ExactHorizonVelcity( obj, x2, obj.outputFile.outputTime(n) );
end
DrawSubplot2d(numel(ind2), obj.outputFile.outputTime, etaNum, etaExt, '$\xi$ (m)');
% DrawSubplot2d(numel(ind2), obj.outputFile.outputTime, uExt, uNum, '$u$ (m/s)');

end

function DrawSubplot2d(Nx2, time, etaNum, etaExt, label)
figure;
for n = 1:Nx2
    subplot( Nx2, 1, n );
    hold on
    plot( time, etaNum(n, :), 'b-', 'LineWidth', 2 );
    plot( time, etaExt(n, :), 'r--', 'LineWidth', 2 );
    box on; grid on;
    xlabel('Time (s)', 'FontSize', 12, 'Interpreter', 'Latex');
    ylabel(label, 'FontSize', 12, 'Interpreter', 'Latex');
end
subplot( Nx2, 1, 1 );
legend({'Numerical results', 'Analytical results'}, ...
    'FontSize', 12, ...
    'Interpreter', 'Latex', ...
    'Location', 'northeastoutside');
end

function eta = ExactWaterLevel( obj, x, time )

Np = numel(x);
Nt = numel(time);
eta = zeros( Np, Nt );
z = zeros(Np, 1);
for n = 1:numel(x)
    temp = obj.matEvaluateExactSolution( x(n), z(n), time );
    eta(n, :) = temp(:);
end

end
