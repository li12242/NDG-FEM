
classdef ConicalLand2d < SWEPreBlanaced2d
    
    properties(Constant)
        %> wet/dry depth threshold
        hmin = 1e-4
        %> gravity acceleration
        gra = 9.8
        %> initial water depth
        h0 = 0.32
        x0 = 12.5
        y0 = 15
        L = 15
        alpha = 0.1
    end
    
    methods
        function obj = ConicalLand2d(N, M, cellType)
            [ mesh ] = makeUniformMesh(N, M, cellType);
%             mesh = readGmshFile( N );
            obj = obj@SWEPreBlanaced2d();
            obj.initPhysFromOptions( mesh );
        end
        
        drawMaxRunup( obj );
        
        function drawGaugeResult( obj )
%             gx = [ 15.71, 15, 15, 15, 15, 15, ...
%                 17.58, 17.88, 18.2, 18.55, 19.63, 21.13, ...
%                 15, 15, 15, 15, 15, 15 ];
%             gy = [ 6.85, 8.4, 9.4, 9.8, 10.12, 10.4, ...
%                 13, 13, 13, 13, 13, 13, ...
%                 15.6, 15.88, 16.2, 16.63, 17.63, 19.13 ];
           
%             gx = [ 6.82, 9.36, 10.36, 12.96, 15.56 ];
%             gy = [ 13.05, 13.80, 13.80, 11.22, 13.80 ];
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
            legend({'Measured data', 'p=1 with coarse mesh', ...
                'p=2 with coarse mesh', 'p=1 with fine mesh'}, ...
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
            
            function drawGaugeValue( measuredFile, time, eta )
                mdata = load( measuredFile );
                plot( mdata(:, 1), mdata(:, 2), 'ro' ); hold on;
                plot( time, eta, 'c' ); 
                xlim([ 5, 20 ]); 
                fontSize = 16;
                xlabel('t (s)', 'Interpreter', 'latex', 'FontSize', fontSize)
                ylabel('$\eta$ (m)', 'Interpreter', 'latex', 'FontSize', fontSize)
            end
        end
        
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
    end
    
    methods(Access=protected)
        function matUpdateExternalField( obj, time, fphys )
            beta = ( obj.h0/obj.L )^2;
            tx = sqrt( 3 * obj.alpha / 4 / beta * ( 1 + obj.alpha ) ) ...
                * sqrt( obj.gra * obj.h0 ) / obj.L * ( time );
%             alpha = obj.alpha;
%             L = obj.L;
%             h0 = obj.h0;
%             xi = sqrt( 3 * alpha * ( 1+alpha ) * L^2 / 4 / h0^2 );
%             tx = xi * sqrt( obj.gra * h0 / L ) * ( time );
            
            sech2 = ( 2 / ( exp(tx) + exp(-tx) ) )^2;
            theta = 1.9;
            h = obj.h0 + obj.alpha * obj.h0 * sech2 * theta;
            for m = 1:obj.Nmesh
                obj.fext{m}(:,:,1) = h;
            end
        end
        
        function fphys = setInitialField( obj )
            fphys = cell( obj.Nmesh, 1 );
            xt = obj.x0;
            yt = obj.y0;
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                fphys{m} = zeros( mesh.cell.Np, mesh.K, obj.Nfield );
                
                r = sqrt( (mesh.x - xt).^2 + (mesh.y - yt).^2 );
                bot = min( 0.625, 0.9 - r/4 );
                bot( bot <= 0 ) = 0.0;
                fphys{m}(:,:,4) = bot;
                temp = obj.h0 - bot;
                temp( temp < 0 ) = 0.0;
                fphys{m}(:,:,1) = temp;
            end
        end
        
        function [ option ] = setOption( obj, option )
            ftime = 20;
            outputIntervalNum = 100;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('temporalDiscreteType') = NdgTemporalIntervalType.DeltaTime;
            option('obcType') = NdgBCType.None;
            option('outputIntervalType') = NdgIOIntervalType.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputNetcdfCaseName') = mfilename;
            option('temporalDiscreteType') = NdgTemporalDiscreteType.RK45;
            option('limiterType') = NdgLimiterType.Vert;
            option('equationType') = NdgDiscreteEquationType.Strong;
            option('integralType') = NdgDiscreteIntegralType.QuadratureFree;
            option('CoriolisType') = CoriolisType.None;
            option('WindType') = WindType.None;
            option('FrictionType') = FrictionType.None;
        end
    end
    
end

function [ mesh ] = readGmshFile( N )
[ path, name, ext ] = fileparts( mfilename('fullpath') );
obj.gmshFile = [ path, '/mesh/conicalLand.msh' ];
mesh = makeGmshFileUMeshUnion2d( N, obj.gmshFile );
end

function [ mesh ] = makeUniformMesh(N, M, type)
bctype = [...
    NdgEdgeType.ZeroGrad, ...
    NdgEdgeType.ZeroGrad, ...
    NdgEdgeType.ClampedDepth, ...
    NdgEdgeType.ZeroGrad];

if (type == NdgCellType.Tri)
    mesh = makeUniformTriMesh(N, [0, 25], [0, 30], M, M, bctype);
elseif(type == NdgCellType.Quad)
    mesh = makeUniformQuadMesh(N, [0, 25], [0, 30], M, M, bctype);
else
    msgID = [mfile, ':inputCellTypeError'];
    msgtext = 'The input cell type should be NdgCellType.Tri or NdgCellType.Quad.';
    ME = MException(msgID, msgtext);
    throw(ME);
end
end% func

