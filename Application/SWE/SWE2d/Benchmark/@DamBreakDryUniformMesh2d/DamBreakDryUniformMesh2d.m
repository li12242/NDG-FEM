classdef DamBreakDryUniformMesh2d < SWEWDPreBlanaced2d
    
    properties( SetAccess = protected )
        theta
    end
    
    properties(Constant)
        %> wet/dry depth threshold
        hmin = 1e-4
        %> gravity acceleration
        gra = 9.8
        %> Dam position
        damPosition = 500
        %> initial water depth
        h0 = 10
    end
    
    methods
        function obj = DamBreakDryUniformMesh2d(N, M, cellType, theta)
            [ mesh ] = makeUniformMesh(N, M, cellType, theta);
            obj = obj@SWEWDPreBlanaced2d();
            obj.theta = theta;
            obj.initPhysFromOptions( mesh );
            obj.fphys = obj.matEvaluatePostFunc( obj.fphys );
        end
        
        function CheckSection( obj )
            Ng = 100;
            xg = linspace(0, 1000, Ng)'; yg = zeros(Ng, 1);
            pos = Analysis2d( obj, xg, yg );
            fInterp = pos.InterpGaugeResult( obj.fphys );
            ftime = obj.getOption('finalTime');
            fext = obj.getExactFunction( ftime );
            fextInterp = pos.InterpGaugeResult( fext );
            
            figure('Color', 'w');
            subplot( 3, 1, 1 ); hold on; grid on; box on;
            plot( xg, fInterp(:,1), 'b.-' );
            plot( xg, fextInterp(:,1), 'r.-' );
            subplot( 3, 1, 2 ); hold on; grid on; box on;
            plot( xg, fInterp(:,2), 'b.-' );
            plot( xg, fextInterp(:,2), 'r.-' );
            subplot( 3, 1, 3 ); hold on; grid on; box on;
            uphyInterp = fInterp(:,2)./fInterp(:,1);
            uextInterp = fextInterp(:,2)./fextInterp(:,1);
            plot( xg, uphyInterp, 'b.-' );
            plot( xg, uextInterp, 'r.-' );
            xlim([0, 1000])
        end
    end
    
    methods(Access=protected)
        function fphys = setInitialField( obj )
            fphys = getExactFunction(obj, 0);
        end
        
        function [ option ] = setOption( obj, option )
            ftime = 20;
            outputIntervalNum = 50;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('outputIntervalType') = enumOutputInterval.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputCaseName') = mfilename;
            option('temporalDiscreteType') = enumTemporalDiscrete.RK45;
            option('limiterType') = enumLimiter.Vert;
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
            option('NumFluxType') = enumSWENumFlux.HLL;
        end
        
        fphys = getExactFunction( obj, time )
    end
    
end

function [ mesh, theta ] = makeUniformMesh(N, M, type, theta)
bctype = [...
    enumBoundaryCondition.ZeroGrad, ...
    enumBoundaryCondition.ZeroGrad, ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall];

if (type == enumStdCell.Tri)
    mesh = makeUniformRotationTriMesh(N, [0, 1000], [-10, 10], ...
        M, ceil(M/50), bctype, 500, 0, theta);
elseif(type == enumStdCell.Quad)
    mesh = makeUniformRotationQuadMesh(N, [0, 1000], [-10, 10], ...
        M, ceil(M/50), bctype, 500, 0, theta );
else
    msgID = [mfile, ':inputCellTypeError'];
    msgtext = 'The input cell type should be NdgCellType.Tri or NdgCellType.Quad.';
    ME = MException(msgID, msgtext);
    throw(ME);
end
end% func

