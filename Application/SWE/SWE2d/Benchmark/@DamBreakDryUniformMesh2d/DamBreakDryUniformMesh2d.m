classdef DamBreakDryUniformMesh2d < SWEConventional2d
    
    properties( SetAccess = protected )
        theta
    end
    
    properties(Constant)
        %> wet/dry depth threshold
        hmin = 1e-5
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
            obj = obj@SWEConventional2d();
            obj.theta = theta;
            obj.initPhysFromOptions( mesh );
            obj.fphys = obj.matEvaluatePostFunc( obj.fphys );
        end
        
        CheckSection( obj, varargin )
        GetNormErr( obj )
        [ fext ] = GetExactFunction( obj, time )
    end
    
    methods(Access=protected)
        function fphys = setInitialField( obj )
            fphys = GetExactFunction(obj, 0);
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
        M, 3, bctype, 500, 0, theta);
elseif(type == enumStdCell.Quad)
    mesh = makeUniformRotationQuadMesh(N, [0, 1000], [-10, 10], ...
        M, 3, bctype, 500, 0, theta );
else
    msgID = [mfile, ':inputCellTypeError'];
    msgtext = 'The input cell type should be NdgCellType.Tri or NdgCellType.Quad.';
    ME = MException(msgID, msgtext);
    throw(ME);
end
end% func

