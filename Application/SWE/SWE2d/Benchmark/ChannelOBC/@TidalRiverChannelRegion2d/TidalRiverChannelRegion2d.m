classdef TidalRiverChannelRegion2d < TidalRiverChannel2d
    
    properties( Constant )
        %> left position of channel
        ChLeft = - 28e3 * 2; % domain length = wave length * lambda
        %> right position of channel
        ChRight = - 28e3 * 1;
    end
    
    properties
        LeftResult
        RightResult
        largeResult
        largeTime
    end
    
    methods( Access = protected, Static  )
        %> set open boundary condition
        function obtype = setOpenBoundaryCondition( )
            obtype = [ NdgEdgeType.Flather, NdgEdgeType.NonLinearFlatherDepth ];
        end
    end
    
    methods( Access = public )
        function obj = TidalRiverChannelRegion2d( N, M, largeSolver )
            obj = obj@TidalRiverChannel2d( N, M );
            
            obj.largeResult = NdgPostProcess( ...
                largeSolver.meshUnion, ...
                'TidalRiverChannelLarge2d.1');
            obj.largeTime = obj.largeResult.time{1};
            obj.AccessLargeResult();
        end
        
        function checkGaugePoint( obj )
            for fld = 1:3
                Ng = 3;
                xg = linspace( obj.ChLeft, obj.ChRight, Ng )';
                yg = obj.ChWidth / 2 * ones( Ng, 1 );

                gaugeLargeResult = ...
                    obj.largeResult.interpolateOutputResultToGaugePoint( xg, yg, yg );

                figure;
                for n = 1:Ng % plot large depth
                    subplot( Ng, 1, n );
                    plot( obj.largeTime, gaugeLargeResult(:, fld, n), 'r.-', ...
                        'LineWidth', 2 );
                    grid on; hold on; box on;
                    xlabel('Time (s)', 'Interpreter', 'latex','FontSize', 16);
                    if fld == 1
                        ylabel('Depth (m)', 'Interpreter', 'latex','FontSize', 16);
                    else
                        ylabel('Flux ($\mathrm{m}^2/\mathrm{s}$)', ...
                            'Interpreter', 'latex','FontSize', 16);
                    end
                end

                pos = makeNdgPostProcessFromNdgPhys( obj );
                gaugeResult = pos.interpolateOutputResultToGaugePoint( xg, yg, yg );
                for n = 1:Ng % plot large depth
                    subplot( Ng, 1, n );
                    plot( pos.time{1}, gaugeResult(:, fld, n), 'b.-', ...
                        'LineWidth', 1 );
                end
                legend({'Large', 'Regional'}, 'box', 'off', ...
                    'Location', 'NorthWest', ...
                    'FontSize', 16, ...
                    'Interpreter', 'latex');
            end
        end
    end
    
    methods( Access = private )
        function AccessLargeResult( obj )
            vbLeftNodeX = obj.meshUnion.x( :, 1 );
            vbLeftNodeY = obj.meshUnion.y( :, 1 );
            obj.LeftResult = ...
                obj.largeResult.interpolateOutputResultToGaugePoint ...
                (vbLeftNodeX, vbLeftNodeY, vbLeftNodeY);
            
            vbRightNodeX = obj.meshUnion.x( :, end );
            vbRightNodeY = obj.meshUnion.y( :, end );
            obj.RightResult = ...
                obj.largeResult.interpolateOutputResultToGaugePoint ...
                (vbRightNodeX, vbRightNodeY, vbRightNodeY);
        end
        
        function mesh = makeUniformMesh( obj, N, M )
            obtype = obj.setOpenBoundaryCondition();
            bctype = [ NdgEdgeType.SlipWall, NdgEdgeType.SlipWall, ...
                obtype(1), obtype(2) ];
            mesh = makeUniformQuadMesh(N, ...
                [-obj.ChLength, 0], [0, obj.ChWidth], M, 1, bctype);
        end% func
    end
    
    methods( Access = protected )
        function fphys = setInitialField( obj )
            fphys = cell( obj.Nmesh, 1 );
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                fphys{m} = zeros( mesh.cell.Np, mesh.K, obj.Nfield );
                fphys{m}(:, :, 1) = + obj.H;
                fphys{m}(:, :, 4) = - obj.H;
            end
        end
        
        function matUpdateExternalField( obj, time, fphys )
            % find previous input boundary time step
            [ idP ] = find( obj.largeTime < time, 1, 'last' );
            % find next input boundary time step
            [ idN ] = find( obj.largeTime > time, 1, 'first' );
            
            if isempty( idP )
                idP = 1;
                delta_time_P = 0;
            else
                delta_time_P = abs( time - obj.largeTime(idP) );
            end
            
            if isempty( idN )
                idN = numel( obj.largeTime );
                delta_time_N = 0;
            else
                delta_time_N = abs( time - obj.largeTime(idN) );
            end
            
            alpha_P = delta_time_P / ( delta_time_P + delta_time_N );
            alpha_N = delta_time_N / ( delta_time_P + delta_time_N );
            
            obj.AccessBoundaryLargeResult( idP, idN, alpha_P, alpha_N );
        end
        
        function [ option ] = setOption( obj, option )
            ftime = 1.3 * ( obj.ChRight - obj.ChLeft ) * 9 / sqrt( obj.gra * obj.H );
            outputIntervalNum = 500;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('obcType') = NdgBCType.None;
            option('outputIntervalType') = NdgIOIntervalType.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputNetcdfCaseName') = ...
                [mfilename, '.', num2str(obj.meshUnion.cell.N)];
            option('temporalDiscreteType') = NdgTemporalDiscreteType.RK22;
            option('limiterType') = NdgLimiterType.None;
            option('equationType') = NdgDiscreteEquationType.Strong;
            option('integralType') = NdgDiscreteIntegralType.GaussQuadrature;
            option('CoriolisType') = SWECoriolisType.None;
            option('WindType') = SWEWindType.None;
            option('FrictionType') = SWEFrictionType.None;
        end
    end
    
    methods( Access = private )
        function AccessBoundaryLargeResult( obj, idP, idN, alphaP, alphaN )
            LeftBoundaryResultP = obj.LeftResult( idP, :, : );
            LeftBoundaryResultN = obj.LeftResult( idN, :, : );
            
            RightBoundaryResultP = obj.RightResult( idP, :, : );
            RightBoundaryResultN = obj.RightResult( idN, :, : );
            
            for m = 1:obj.Nmesh
                for fld = 1:3
                    temp = alphaP .* LeftBoundaryResultP(:, fld, :) + ...
                        alphaN .* LeftBoundaryResultN(:, fld, :);
                    obj.fext{m}(1, 1, fld) = temp(:, :, 1);
                    obj.fext{m}(2, 1, fld) = temp(:, :, 2);
                    obj.fext{m}(3, 1, fld) = temp(:, :, 3);
                    obj.fext{m}(4, 1, fld) = temp(:, :, 4);
                    
                    
                    temp = alphaP .* RightBoundaryResultP(:, fld, :) + ...
                        alphaN .* RightBoundaryResultN(:, fld, :);
                    obj.fext{m}(1, end, fld) = temp(:, :, 1);
                    obj.fext{m}(2, end, fld) = temp(:, :, 2);
                    obj.fext{m}(3, end, fld) = temp(:, :, 3);
                    obj.fext{m}(4, end, fld) = temp(:, :, 4);
                end
            end
        end
    end
    
end

