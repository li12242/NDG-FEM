function AnalysisResult3d( obj )
%ANALYSISRESULT3D Summary of this function goes here
%   Detailed explanation goes here

OutstepNum = obj.outputFile.outputStep;
Nz = obj.mesh3d.Nz;
Npz = obj.mesh3d.cell.Npz;
Nph = obj.mesh3d.cell.Nph;
Nptol = obj.mesh3d.cell.Np * obj.mesh3d.K;

% analysis vertical velocity
cellId2 = 21;
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
z = obj.mesh3d.z( ind ) * obj.H;
Nnode = numel( z );

wNum = zeros( Nnode, OutstepNum );
wExt = zeros( Nnode, OutstepNum );

w_ind = ind + 2 * Nptol;
for n = 1:OutstepNum
    [ ~, fphys3d ] = obj.outputFile.readOutputResult(n);
    
    wNum( :, n ) = fphys3d( w_ind );
    wExt( :, n ) = ExactVerticalVelcity( obj, x, z, obj.outputFile.outputTime(n) );
end

% DrawSubplot2d( Nnode, obj.outputFile.outputTime, wNum, wExt, '$w$ (m/s)' );

T = 2 * pi / obj.omega;
Tstep = T ./ obj.outputTimeInterval;
Nt = 4;
timeStep = round( linspace( 1, Tstep, Nt ) );
DrawSubplot3d( Nt, timeStep, z, wNum, wExt );

end

function DrawVerticalAnimation( Nt, wNum, wExt )
figure;
hold on
p1 = plot( wNum(:, t), z, 'b-', 'LineWidth', 2 );
p2 = plot( wExt(:, t), z, 'r--', 'LineWidth', 2 );
box on; grid on;
xlabel('$\omega$ (m/s)', 'FontSize', 12, 'Interpreter', 'Latex');
ylabel('$z$ (m)', 'FontSize', 12, 'Interpreter', 'Latex');
ylim([-12, 0]);

for n = 1:Nt
    
end
end

function DrawSubplot2d(Nx2, time, wNum, wExt, label)
figure
for n = 1:Nx2
    subplot( Nx2, 1, n );
    hold on
    plot( time, wNum(n, :), 'b-', 'LineWidth', 2 );
    plot( time, wExt(n, :), 'r--', 'LineWidth', 2 );
    box on; grid on;
    ylim([ -2.5, 2.5 ] * 1e-4);
    xlabel('Time (s)', 'FontSize', 12, 'Interpreter', 'Latex');
    ylabel(label, 'FontSize', 12, 'Interpreter', 'Latex');
end
subplot( Nx2, 1, 1 );
legend({'Numerical results', 'Analytical results'}, ...
    'FontSize', 12, ...
    'Interpreter', 'Latex', ...
    'Location', 'northeastoutside');
end

function DrawSubplot3d( Nt, timeStep, z, wNum, wExt )
figure;
hold on
for n = 1 : Nt
    t = timeStep(n);
    plot( wNum(:, t), z, 'b-o', 'LineWidth', 2, ...
        'MarkerSize', 5, 'MarkerFaceColor', 'b' );
    plot( wExt(:, t), z, 'r--', 'LineWidth', 2 );
    box on; grid on;
    xlabel('$\omega$ (m/s)', 'FontSize', 12, 'Interpreter', 'Latex');
    ylabel('Depth (m)', 'FontSize', 12, 'Interpreter', 'Latex');
    ylim([-12, 0]);
    %xlim([ -4, 4 ] * 1e-5);
end
legend({'Numerical', 'Analytical'}, ...
    'FontSize', 12, ...
    'Interpreter', 'Latex', ...
    'Location', 'northeast');
end

function w = ExactVerticalVelcity( obj, x, z, time )

sgH = sqrt( obj.gra * obj.H );
w = - obj.a * obj.omega * (z + obj.H) / obj.H / ...
    cos( obj.omega * obj.ChLength / sgH ) .* ...
    cos( obj.omega * ( obj.ChLength - x ) / sgH ) * ...
    sin( obj.omega * time );

end
