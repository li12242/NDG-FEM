classdef SteadyMount2d < SWEPreBlanaced2d & CSBAbstractTest
     
    properties( Constant )
        %> wet/dry depth threshold
        hmin = 1e-2
        %> gravity acceleration
        gra = 9.8
        
        x0 = 0.5
        y0 = 0.5
    end
    
    properties
        N
    end
    
    methods
        function obj = SteadyMount2d( N, M, cellType )
            [ mesh ] = makeUniformMesh(N, M, cellType);
            obj = obj@SWEPreBlanaced2d();
            obj.N = N;
            obj.initPhysFromOptions( mesh );
            obj.fphys = obj.matEvaluatePostFunc( obj.fphys );
        end
        
        function drawFluxError( obj )
            ncfile = [obj.getOption('outputNetcdfCaseName'), '.1-1.nc'];
            time = ncread( ncfile, 'time' );
            pos = makeNdgPostProcessFromNdgPhys( obj );
            err = zeros( pos.Nt, 1 );
            for t = 1:numel( time )
                fext = obj.getExactFunction( time(t) );
                fphys = pos.accessOutputResultAtStepNum( t );
                temp = pos.evaluateNormErr2( fphys, fext );
                err( t,1 ) = temp( 2 );
                %err( t,2 ) = temp( 3 );
            end
            hold on;  grid on;
            plot( time, err(:,1), 'g', 'LineWidth', 1.5 );
            xlabel('$t$ (s)', 'FontSize', 16, 'Interpreter', 'Latex');
            ylabel('$\left\| hu \right\|_2$', 'FontSize', 16, 'Interpreter', 'Latex');
            legend( {['$p=', num2str(obj.N), '$']}, ...
                'FontSize', 16, 'Interpreter', 'Latex', 'box', 'off',...
                'Location', 'NorthEast')
        end
        
        function coutourFlux( obj )
            pos = makeNdgPostProcessFromNdgPhys( obj );
            fphys = pos.accessOutputResultAtStepNum( pos.Nt );
            Ng = 50;
            tol = 0;
            xg = linspace( 0+tol, 1-tol, Ng );
            yg = linspace( 0+tol, 1-tol, Ng );
            [ x, y ] = meshgrid( xg, yg );
            [ fg ] = pos.interpolatePhysFieldToGaugePoint(...
                fphys, x(:), y(:), y(:) );
            flux = reshape( fg(:, 2), Ng, Ng );
            contourf( xg, yg, real(flux), 20 );
            colorbar;
            colormap jet;
%             set( gca, 'CLim', [-1.5, 1.5]*1e-4 );
            xlabel('$x$(m)', 'FontSize', 16, 'Interpreter', 'Latex');
            ylabel('$y$(m)', 'FontSize', 16, 'Interpreter', 'Latex');
        end
            
    end
    methods(Access=protected)
        
        function fphys = setInitialField( obj )
            fphys = getExactFunction(obj, 0);
        end
        
        function [ option ] = setOption( obj, option )
            ftime = 10;
            outputIntervalNum = 50;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('outputIntervalType') = enumOutputInterval.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputCaseName') = [mfilename, '_', num2str(obj.N)];
            option('temporalDiscreteType') = enumTemporalDiscrete.RK45;
            option('limiterType') = enumLimiter.Vert;
            option('SWELimiterType') = enumSWELimiter.OnElevation;
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
            option('CoriolisType')= enumSWECoriolis.None;
            option('WindType') = enumSWEWind.None;
            option('FrictionType') = enumSWEFriction.None;
        end
        
        function fphys = getExactFunction( obj, time )
            
            fphys = cell( obj.Nmesh, 1 );
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                bot = max( 0, 0.25 - 2.5*(mesh.x - obj.x0).^2 - 2.5*(mesh.y - obj.y0).^2 );
                bmin = min(bot);
                bmean = mesh.GetMeshAverageValue( bot );
                bmean( bmean <= 0 ) = 0;
                theta = bsxfun(@min, 1, bmean./(bmean - bmin) );
                bot = bsxfun(@plus, bsxfun(@times, bsxfun(@plus, bot, -bmean ), theta ), bmean);
% %                 sigma = 25*1e3/(33*33);
% %                 temp = sigma*( (mesh.x - obj.x0).^2 + (mesh.y - obj.y0).^2 );
% %                 bot = exp( -temp );
% %                 fphys{m}(:,:,4) = bot;
% %             end
% %             fphys = limiter.matLimit( fphys, 4 );
% %             for m = 1:obj.Nmesh
% %                 bot = fphys{m}(:,:,4);
                eta = 0.2 * ones( mesh.cell.Np, mesh.K );
                h = max( eta - bot, 0 );
                ind = ( max(h) > 0 );
                bot(:, ind) = eta(:, ind) - h(:, ind);
                fphys{m}(:,:,1) = h;
                fphys{m}(:,:,4) = bot;
            end
        end
    end
    
end

function [ mesh ] = makeUniformMesh(N, M, type)
bctype = [...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall];

if (type == enumStdCell.Tri)
    mesh = makeUniformTriMesh(N, [0, 1], [0, 1], M, M, bctype);
elseif(type == enumStdCell.Quad)
    mesh = makeUniformQuadMesh(N, [0, 1], [0, 1], M, M, bctype);
else
    msgID = [mfile, ':inputCellTypeError'];
    msgtext = 'The input cell type should be NdgCellType.Tri or NdgCellType.Quad.';
    ME = MException(msgID, msgtext);
    throw(ME);
end
end% func
