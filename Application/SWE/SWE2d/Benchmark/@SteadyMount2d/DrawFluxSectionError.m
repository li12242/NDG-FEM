function DrawFluxSectionError( obj )
%DRAWFLUXSECTIONERROR Summary of this function goes here
%   Detailed explanation goes here

figure(1);
LineStyle = 'r-.';
% LineStyle = 'b-';

kid = (obj.M^2 - obj.M + 1):obj.M ^ 2;
fid = 2;
nid = obj.meshUnion.cell.Fmask(:, fid);
% Ntol = obj.meshUnion.K * obj.meshUnion.cell.Np;
xg = obj.meshUnion.x( nid, kid );
obj.matEvaluateRHS( obj.fphys );
hold on;
plot( xg, obj.frhs{1}( nid, kid, 2 ), LineStyle, 'LineWidth', 2 );

xlabel('$x$ (m)', 'FontSize', 16, 'Interpreter', 'Latex');
ylabel('$\frac{\partial }{\partial x} F_h + S_h$', ...
    'FontSize', 16, 'Interpreter', 'Latex');

end

% [legh,objh,outh,outm] = legend({'Conventional', 'Our Scheme'}, ...
%     'Interpreter', 'Latex', 'FontSize', 14);
