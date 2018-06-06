classdef SolitaryRossbyWave < SWEConventional2d
    
    
    properties( SetAccess = protected )
    end
    
    properties(Constant)
        %> wet/dry depth threshold
        hmin = 1e-4
        %> gravity acceleration
        gra = 1
        %> amplitude of the solitary wave
        a = 0.395
        
    end
    
    methods
        function obj = SolitaryRossbyWave(N, M, cellType )
            [ mesh ] = makeUniformMesh(N, M, cellType );
            obj = obj@SWEConventional2d();
            obj.initPhysFromOptions( mesh );
        end
        
    end%methods
    
    methods(Access=protected)
        function fphys = setInitialField( obj )
            fphys = getInitialFunction(obj);
        end
        
        function [ option ] = setOption( obj, option )
            ftime = 40;
            outputIntervalNum = 50;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('temporalDiscreteType') = NdgTemporalIntervalType.DeltaTime;
            option('obcType') = NdgBCType.None;
            option('outputIntervalType') = NdgIOIntervalType.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputNetcdfCaseName') = 'SolitaryRossbyWave';
            option('temporalDiscreteType') = NdgTemporalDiscreteType.RK45;
            option('limiterType') = NdgLimiterType.None;
            option('equationType') = NdgDiscreteEquationType.Strong;
            option('integralType') = NdgDiscreteIntegralType.QuadratureFree;
            option('CoriolisType')=SWECoriolisType.Beta;
            option('CoriolisParameter_f0')=0;
            option('CoriolisParameter_beta')=1;
            option('WindType')=SWEWindType.None;
            option('FrictionType')=SWEFrictionType.None;
        end
        
    end%methods
    
end%classdef

function [ mesh, theta ] = makeUniformMesh(N, M, type, theta)
bctype = [...
    NdgEdgeType.SlipWall, ...
    NdgEdgeType.SlipWall, ...
    NdgEdgeType.SlipWall, ...
    NdgEdgeType.SlipWall];
%SNWE

if (type == NdgCellType.Tri)
    mesh = makeUniformTriMesh(N, [-24, 24], [-8, 8], ...
        M, M, bctype );
    
elseif(type == NdgCellType.Quad)
    
    mesh = makeUniformQuadMesh(N, [-24, 24], [-8, 8], ...
        M, M, bctype );
else
    msgID = 'SolitaryRossbyWave:inputCellTypeError';
    msgtext = 'The input cell type should be NdgCellType.Tri or NdgCellType.Quad.';
    ME = MException(msgID, msgtext);
    throw(ME);
end
end% func