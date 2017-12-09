classdef FlatBathNonlinearStommel < SWEAbstractDBN82d

    
    properties( SetAccess = protected )
    end
    
    properties(Constant)
         %> wet/dry depth threshold
        hmin = 1e-4                
        %> gravity acceleration        
        gra = 10
        
        %> constant part of coriolis parameter
%         f0 = 1e-4
        %> meridional vatiation of coriolis parameter
%         beta = 1e-11
        
        %> bottom friction coefficient
%         r = 2e-6
        
        %> wind stress coefficient
%         cd = 2.6e-3
        %> dencity of sea water(kg/m^3)
%         rou = 1000
        %> dencity of air(kg/m^3)
%         rouair = 1.226
        
    end
    
    methods
        function obj = FlatBathNonlinearStommel(N, M, cellType)
            [ mesh ] = makeUniformMesh(N, M, cellType );
            obj = obj@SWEAbstractDBN82d();
            obj.initPhysFromOptions( mesh );
%             obj.fext = obj.setFinalField();
        end 
        
%         function  setFinalField( obj )
%             ftime = obj.getOption('finalTime');
%             obj.fext = getExactFunction(obj, ftime);%fext的终止场设为精确解的终止场
%         end
        
    end%methods
    
    methods(Access=protected)
        function fphys = setInitialField( obj )
            fphys = getInitialFunction(obj);
        end
        
        function [ option ] = setOption( obj, option )
            ftime = 7430400;
            outputIntervalNum = 50;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('temporalDiscreteType') = NdgTemporalIntervalType.DeltaTime;
            option('obcType') = NdgBCType.None;
            option('outputIntervalType') = NdgIOIntervalType.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputNetcdfCaseName') = 'FlatBathNonlinearStommelProblem';
            option('temporalDiscreteType') = NdgTemporalDiscreteType.RK45;
            option('limiterType') = NdgLimiterType.None;
            option('equationType') = NdgDiscreteEquationType.Strong;
            option('integralType') = NdgDiscreteIntegralType.QuadratureFree;
            option('CoriolisType')=CoriolisType.Beta;
            option('CoriolisParameter_f0')=1e-4;
            option('CoriolisParameter_beta')=1e-11;
            option('WindType')=WindType.Stress;
            option('WindSterssCoefficient_cd')=2.6e-3;
            option('DensityofAir')=1.226;
            option('DensityofWater')=1000;
            option('FrictionType')=FrictionType.Linear;
            option('FrictionCoefficient_r')=2e-6;
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
    mesh = makeUniformTriMesh(N, [0, 1000000], [0, 1000000], ...
        M, M, bctype );

elseif(type == NdgCellType.Quad)

    mesh = makeUniformQuadMesh(N, [0, 1000000], [0, 1000000], ...
        M, M, bctype );
else
    msgID = 'FlatBathNonlinearStommel:inputCellTypeError';
    msgtext = 'The input cell type should be NdgCellType.Tri or NdgCellType.Quad.';
    ME = MException(msgID, msgtext);
    throw(ME);
end
end% func

