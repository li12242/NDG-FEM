function drawResult( obj )
          
            timeFrac = [0.3, 0.4, 0.45, 0.52];
            [ conicalPos ] = makeNdgPostProcessFromNdgPhys( obj );
            [ time ] = ncread( conicalPos.outputFile{1}, 'time' );
            Nt = numel( timeFrac );
            timeStep = zeros(Nt, 1);
            for n = 1:Nt
                timeSpecific = timeFrac(n) * max(time);
                [ ~, timeStep(n) ] = min( abs(time - timeSpecific) );
            end
            % draw 3d plot
            for n = 1:Nt
                figure;
                ts = timeStep(n);
                field = conicalPos.accessOutputResultAtStepNum( ts );
                
                % draw 3d plot
                draw3dBottom( obj.meshUnion, obj.fphys{1}(:,:,4) );
                eta = field{1}(:,:,1) + obj.fphys{1}(:,:,4);
                eta( field{1}(:,:,1) < 1.5e-3 ) = nan;
                draw3dSurface( obj.meshUnion, eta );
                colormap jet;
                
                xlim([12.5 - 10, 12.5 + 10]);
                ylim([15 - 10, 15 + 10]);
                zlim([0.32 - 0.05, 0.32 + 0.1]);
                if n < 3
                    view([-40, 23]);
                else
                    view([77, 23]);
                end
                set( gca, 'CLim', [0.32 - 0.02, 0.32 + 0.06] );
                xlabel('$x$ (m)', 'FontSize', 16, 'Interpreter', 'Latex');
                ylabel('$y$ (m)', 'FontSize', 16, 'Interpreter', 'Latex');
                zlabel('$h - h_0$ (m)', 'FontSize', 16, 'Interpreter', 'Latex');
                colorbar;
            end
            
                function handle = draw3dSurface( mesh, zvar )
                    EToV = ones(mesh.K, 1)*mesh.cell.Fmask(:)';
                    EToV = EToV + ( mesh.cell.Np*(0:mesh.K-1) )' * ones(1, mesh.cell.TNfp);
                    handle = patch(...
                        'Vertices', [mesh.x(:), mesh.y(:), zvar(:)], ...
                        'Faces', EToV, ...
                        'FaceColor', 'interp', ...
                        'EdgeColor', 'none', ...
                        'FaceVertexCData', zvar(:));
                    box on;
                    grid on;
                end

                function handle = draw3dBottom( mesh, zvar )
                    EToV = ones(mesh.K, 1)*mesh.cell.Fmask(:)';
                    EToV = EToV + ( mesh.cell.Np*(0:mesh.K-1) )' * ones(1, mesh.cell.TNfp);
                    handle = patch(...
                        'Vertices', [mesh.x(:), mesh.y(:), zvar(:)], ...
                        'Faces', EToV, ...
                        'FaceColor', [.67, .67, .67], ...
                        'EdgeColor', 'none', ...
                        'FaceVertexCData', zvar(:));
                    box on;
                    grid on;
                end
        end