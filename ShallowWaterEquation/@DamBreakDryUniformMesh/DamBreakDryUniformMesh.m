classdef DamBreakDryUniformMesh < SWESolidTopography
    
    properties( SetAccess = protected )
        theta
    end
    
    properties(Constant)
        %> wet/dry depth threshold
        hmin = 1e-4
        %> Dam position
        damPosition = 500
        %> initial water depth
        h0 = 10
    end
    
    methods
        function obj = DamBreakDryUniformMesh(N, M, cellType, theta)
            [ mesh ] = makeUniformMesh(N, M, cellType, theta);
            obj = obj@SWESolidTopography();
            obj.theta = theta;
            obj.setNdgPhys( mesh );
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
            option('LimiterType') = NdgLimiterType.VertexLimiter;
            option('temporalDiscreteType') = NdgIntervalType.DeltaTime;
            option('obcType') = NdgBCType.None;
            option('outputIntervalType') = NdgIntervalType.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputNetcdfCaseName') = 'DamBreakDryUniformMesh';
            option('temporalDiscreteType') = NdgTemporalDiscreteType.RK45;
        end
        
        fphys = getExactFunction( obj, time )
    end
    
end

function [ mesh, theta ] = makeUniformMesh(N, M, type, theta)
bctype = [...
    NdgEdgeType.ZeroGrad, ...
    NdgEdgeType.ZeroGrad, ...
    NdgEdgeType.SlipWall, ...
    NdgEdgeType.SlipWall];

if (type == NdgCellType.Tri)
    mesh = makeUniformRotationTriMesh(N, [0, 1000], [-10, 10], ...
        M, ceil(M/50), bctype, 500, 0, theta);
elseif(type == NdgCellType.Quad)
    mesh = makeUniformRotationQuadMesh(N, [0, 1000], [-10, 10], ...
        M, ceil(M/50), bctype, 500, 0, theta );
else
    msgID = 'DamBreakDryUniformMesh:inputCellTypeError';
    msgtext = 'The input cell type should be NdgCellType.Tri or NdgCellType.Quad.';
    ME = MException(msgID, msgtext);
    throw(ME);
end
end% func

