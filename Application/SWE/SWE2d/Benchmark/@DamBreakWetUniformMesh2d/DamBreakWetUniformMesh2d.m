classdef DamBreakWetUniformMesh2d < SWEConventional2d
    
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
        %> water depth at the upstream (left)
        h0 = 10
        %> water depth at the downstream (right)
        h1 = 2
    end
    
    methods
        function obj = DamBreakWetUniformMesh2d(N, M, cellType, theta)
            [ mesh ] = makeUniformMesh(N, M, cellType, theta);
            obj = obj@SWEConventional2d();
            obj.theta = theta;
            obj.initPhysFromOptions( mesh );
        end
        
        function CheckSection( obj )
            Ng = 100;
            xg = linspace(0, 1000, Ng)'; yg = zeros(Ng, 1);
            pos = Analysis2d( obj, xg, yg );
            fphyInterp = pos.InterpGaugeResult( obj.fphys );
            ftime = obj.getOption('finalTime');
            fext = obj.getExactFunction( ftime );
            fextInterp = pos.InterpGaugeResult( fext );
            
            figure('Color', 'w'); hold on; grid on; box on;
            plot( xg, fphyInterp(:,1), 'b.-' );
            plot( xg, fextInterp(:,1), 'r.-' );
            figure('Color', 'w'); hold on; grid on; box on;
            plot( xg, fphyInterp(:,2), 'b.-' );
            plot( xg, fextInterp(:,2), 'r.-' );
            figure('Color', 'w'); hold on; grid on; box on;
            uphyInterp = fphyInterp(:,2)./fphyInterp(:,1);
            uextInterp = fextInterp(:,2)./fextInterp(:,1);
            plot( xg, uphyInterp, 'b.-' );
            plot( xg, uextInterp, 'r.-' );
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
            option('temporalDiscreteType') = enumTemproalInterval.DeltaTime;
            option('outputIntervalType') = enumOutputInterval.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputCaseName') = mfilename;
            option('temporalDiscreteType') = enumTemporalDiscrete.RK45;
            option('limiterType') = enumLimiter.Vert;
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
            option('NumFluxType') = enumSWENumFlux.LF;
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
    msgID = 'DamBreakDryUniformMesh:inputCellTypeError';
    msgtext = 'The input cell type should be NdgCellType.Tri or NdgCellType.Quad.';
    ME = MException(msgID, msgtext);
    throw(ME);
end
end% func
