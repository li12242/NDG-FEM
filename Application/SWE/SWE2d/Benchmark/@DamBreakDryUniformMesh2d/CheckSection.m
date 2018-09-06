function CheckSection(obj, varargin)

Ng = 300;
xg = linspace(0, 1000, Ng)'; 
yg = zeros(Ng, 1);

% time = [4, 8, 12];
time = [4, 12, 20];
% dt = 0.15;
dt = 0;
pos = Analysis2d( obj, xg(:), yg(:) );

marker = {'s', 'o', '^'};
color = {[1, 0, 0], [0, 0, 1], [0.2, 0.9, 0.2]};
for n = 1:numel(time)
    step = find( obj.outputFile.outputTime >= time(n), 1 );
    fphys{1} = obj.outputFile.readOutputResult( step );
    fext = obj.GetExactFunction( obj.outputFile.outputTime(step) - dt );
    
    fphyInterp = pos.InterpGaugeResult( fphys );
    fextInterp = pos.InterpGaugeResult( fext );

    figure(1); hold on;
    plot( xg, fextInterp(:,1), 'k',  'LineWidth', 1.5, ...
        'DisplayName', ['analytical h ', num2str(n)] ); 
    plot( xg, fphyInterp(:,1), marker{n}, 'Color', color{n}, 'LineWidth', 2, ...
        'DisplayName', ['numerical h ', num2str(n)] ); 
    box on; grid on;
    xlabel('$x$ (m)', 'FontSize', 16, 'Interpreter', 'Latex');
    ylabel('Depth (m)', 'FontSize', 16, 'Interpreter', 'Latex');
    if( nargin > 1 )
        xlim([500, 800]); ylim([0, 2.5])
    else
        xlim([0, 1000]); 
    end

    figure(2); hold on;
    plot( xg, fextInterp(:,2), 'k', 'LineWidth', 1.5, ...
        'DisplayName', ['analytical hu ', num2str(n)] ); 
    plot( xg, fphyInterp(:,2), marker{n}, 'Color', color{n}, 'LineWidth', 2, ...
        'DisplayName', ['numerical hu ', num2str(n)] ); 
    xlabel('$x$ (m)', 'FontSize', 16, 'Interpreter', 'Latex');
    ylabel('Flux ($\mathrm{m}^2/$s)', 'FontSize', 16, 'Interpreter', 'Latex');
    box on; grid on;
    if( nargin > 1 )
        xlim([500, 800]); ylim([0, 25]);
    else
        xlim([0, 1000]); ylim([0, 35]);
    end
    
    figure(3); hold on;
    uphyInterp = fphyInterp(:,2)./fphyInterp(:,1);
    uextInterp = fextInterp(:,2)./fextInterp(:,1);
    plot( xg, uextInterp, 'k', 'LineWidth', 1.5);
    plot( xg, uphyInterp, marker{n}, 'Color', color{n} ); 
    xlim([0, 1000]); box on; grid on;
end

AddLegend(1)
AddLegend(2)

end

function AddLegend(n)
figure(n);
[~, t2] = legend({'Exact', '$t = 4$ s', '$t = 12$ s', '$t = 20$ s'}, ...
'FontSize', 16, 'Interpreter', 'Latex', 'box', 'off');
for n = 1:4
t2(n).FontSize = 16;
t2(n).Interpreter = 'Latex';
end

t2(9).Visible = 'off';
t2(10).Marker = 'o';
t2(10).Color = 'b';
t2(10).LineWidth = 0.5;
t2(12).Marker = '^';
t2(12).Color = [0.2, 0.9, 0.2];
end