function drawGaugeResult( obj )
            gx = [ 6.36, 8.9, 9.9, 12.5, 15.1 ];
            gy = [ 14.25, 15, 15, 12.42, 15 ];
            timeDelta = 2.6;
            [ conicalPos ] = makeNdgPostProcessFromNdgPhys( obj );
            [ result ] = conicalPos.interpolateOutputResultToGaugePoint( gx, gy, gx );
            [ time ] = ncread( conicalPos.outputFile{1}, 'time' ) + timeDelta;
            [ gaugeValue ] = conicalPos.interpolatePhysFieldToGaugePoint( ...
                obj.fphys, gx, gy, gx );
            [ bot ] = gaugeValue(:, 4)';
            [ path, name, ext ] = fileparts( mfilename('fullpath') );
            
            % gauge point #3
            figure(1); set(gcf, 'Position', [440, 630, 560, 160]);
            ind = 1;
            temp = result(ind, 1, :) - obj.h0 + bot(ind);
            drawGaugeValue( [path, '/result/G31.csv'], time, temp(:) );
            legend({'Measured data', 'p=1 RKDG', 'p=2 RKDG'}, ...
                    'Interpreter', 'latex', ...
                    'FontSize', 16, 'box', 'off');
            % gauge point #6
            figure(2); set(gcf, 'Position', [440, 630, 560, 160]);
            ind = 2;
            temp = result(ind, 1, :) - obj.h0 + bot(ind);
            drawGaugeValue( [path, '/result/G61.csv'], time, temp(:) );

            % gauge point #9
            figure(3); set(gcf, 'Position', [440, 630, 560, 160]);
            ind = 3;
            temp = result(ind, 1, :) - obj.h0 + bot(ind);
            drawGaugeValue( [path, '/result/G91.csv'], time, temp(:) );
            
            % gauge point #16
            figure(4); set(gcf, 'Position', [440, 630, 560, 160]);
            ind = 4;
            temp = result(ind, 1, :) - obj.h0 + bot(ind);
            drawGaugeValue( [path, '/result/G161.csv'], time, temp(:) );
            
            % gauge point #22
            figure(5); set(gcf, 'Position', [440, 630, 560, 160]);
            ind = 5;
            temp = result(ind, 1, :) - obj.h0 + bot(ind);
            drawGaugeValue( [path, '/result/G221.csv'], time, temp(:) );
            
            
        end

        function drawGaugeValue( measuredFile, time, eta )
                mdata = load( measuredFile );
                plot( mdata(:, 1), mdata(:, 2), 'ro' ); hold on;
                plot( time, eta, 'c' ); 
                xlim([ 5, 20 ]); 
                fontSize = 16;
                xlabel('t (s)', 'Interpreter', 'latex', 'FontSize', fontSize)
                ylabel('$\eta$ (m)', 'Interpreter', 'latex', 'FontSize', fontSize)
            end