function AnalysisResult3d( obj )
%ANALYSISRESULT3D Summary of this function goes here
%   Detailed explanation goes here

OutstepNum = obj.outputFile.outputStep;
Nz = obj.mesh3d.Nz;
Npz = obj.mesh3d.cell.Npz;
Nph = obj.mesh3d.cell.Nph;

% analysis vertical velocity
cellId2 = 11;
cellId3 = (cellId2 - 1) * Nz + (Nz:-1:1);
%cellId3 = (cellId2 - 1) * Nz + 1;
nodeId3 = ((1:Npz) - 1) * Nph + 1;
sk = 1;
for n = 1:numel( cellId3 )
    for m = 1:numel( nodeId3 )
        ind(sk) = sub2ind( [obj.mesh3d.cell.Np, obj.mesh3d.K], ...
            nodeId3(m), cellId3(n) );
        sk = sk + 1;
    end
end
x = obj.mesh3d.x( ind );
sigma = obj.mesh3d.z( ind );
Nnode = numel( sigma );
uNum = zeros( Nnode, OutstepNum );

Nptol = obj.mesh3d.cell.Np * obj.mesh3d.K;
u_ind = ind;
for n = 1:OutstepNum
    [ ~, fphys3d ] = obj.outputFile.readOutputResult(n);
    uNum( :, n ) = fphys3d( u_ind );
end

% DrawSubplot2d( Nnode, obj.outputFile.outputTime, wNum, wExt, '$w$ (m/s)' );

Tstep = ceil( 2*pi/obj.omega/obj.outputTimeInterval );
Nt = 5;
timeStep = round( linspace( 1.2 * Tstep, OutstepNum, Nt ) );
uExt = ExactVerticalVelcity( obj, x, sigma, obj.outputFile.outputTime(timeStep) );

z = sigma .* obj.fphys2d{1}(1, cellId2, 4); % change to z coordinate
DrawSubplot3d( Nt, z, uNum(:, timeStep), uExt );

end

function DrawSubplot3d( Nt, z, uNum, uExt )
figure;
hold on
for n = 1 : Nt
    plot( uNum(:, n), z, 'b-o', 'LineWidth', 2, 'MarkerSize', 5 );
    plot( uExt(:, n), z, 'r--', 'LineWidth', 2 );
    box on; grid on;
    xlabel('$\omega$ (m/s)', 'FontSize', 12, 'Interpreter', 'Latex');
    ylabel('Depth (m)', 'FontSize', 12, 'Interpreter', 'Latex');
    %ylim([-12, 0]);
    %xlim([ -4, 4 ] * 1e-5);
end
ylim([-22.5, 0]);
legend({'$\tau = 10^{-3}$', 'Analytical'}, ...
    'FontSize', 12, ...
    'Interpreter', 'Latex', ...
    'Location', 'Best');
end

function u = ExactVerticalVelcity( obj, x, z, time )

Np = numel(x);
Nt = numel(time);
u = zeros( Np, Nt );
for n = 1:numel(x)
    [~, temp] = obj.matEvaluateExactSolution( x(n), z(n), time );
    u(n, :) = temp(:);
end

end
