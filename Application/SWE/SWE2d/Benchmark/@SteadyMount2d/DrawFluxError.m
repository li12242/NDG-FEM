function DrawFluxError( obj )
time = obj.outputFile.outputTime;
err = zeros( obj.outputFile.outputStep, 1 );
for t = 1 : obj.outputFile.outputStep
    fext = obj.getExactFunction( time(t) );
    fphys = obj.outputFile.readOutputResult( t );
    temp = obj.meshUnion.evaluateNormErr2( fphys, fext{1} );
    err( t ) = temp( 2 ) + temp(3);
end
hold on;  grid on; box on;

err = err * 1e-9;
% plot( time(1:(end-1) ), 0.15*(err(1:(end-1) ) - err_m) + 2.223e-17, 'r-.', 'LineWidth', 1.5 );
plot( time(1:(end-1) ), err(1:(end-1) ), 'b', 'LineWidth', 1.5 );
xlabel('$t$ (s)', 'FontSize', 16, 'Interpreter', 'Latex');
ylabel('$\left\| hu \right\|_2$', 'FontSize', 16, 'Interpreter', 'Latex');
% legend( {['$p=', num2str(obj.N), '$']}, ...
%     'FontSize', 16, 'Interpreter', 'Latex', 'box', 'off',...
%     'Location', 'NorthEast')
end

