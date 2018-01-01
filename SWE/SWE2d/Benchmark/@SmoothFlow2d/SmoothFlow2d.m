classdef SmoothFlow2d < SWEAbstractDBN62d
    %SMOOTHFLOW2D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
        hmin = 1e-4;
        gra = 9.81;
    end
    properties
        
    end
    
    methods
        function obj = SmoothFlow2d(N, M, cellType)
            [ mesh ] = makeUniformMesh(N, M, cellType);
            obj = obj@SWEAbstractDBN62d();
            obj.initPhysFromOptions( mesh );
        end
    end
    
    methods(Access=protected)
        function fphys = setInitialField( obj )
            fphys = cell( obj.Nmesh, 1 );
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                fphys{m} = zeros( mesh.cell.Np, mesh.K, obj.Nfield );
                fphys{m}(:,:,1) = 10 + exp( sin( 2*pi*mesh.x ) ).*cos(2*pi*mesh.y);
                fphys{m}(:,:,2) = sin( cos(2*pi*mesh.x) ).*sin(2*pi*mesh.y);
                fphys{m}(:,:,3) = cos( 2*pi*mesh.x ).*cos(sin(2*pi*mesh.y));
                fphys{m}(:,:,4) = sin( 2*pi*mesh.x ) + cos(2*pi*mesh.y);
            end
        end
        
        function [ option ] = setOption( obj, option )
            ftime = 0.02;
            outputIntervalNum = 10;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('temporalDiscreteType') = NdgTemporalIntervalType.Constant;
            option('obcType') = NdgBCType.None;
            option('outputIntervalType') = NdgIOIntervalType.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputNetcdfCaseName') = mfilename;
            option('temporalDiscreteType') = NdgTemporalDiscreteType.RK45;
            option('limiterType') = NdgLimiterType.None;
            option('equationType') = NdgDiscreteEquationType.Strong;
            option('integralType') = NdgDiscreteIntegralType.QuadratureFree;
            option('CoriolisType')= CoriolisType.None;
            option('WindType') = WindType.None;
            option('FrictionType') = FrictionType.None;
            option('WellBlancedType') = true;
        end
    end
    
end

function [ mesh ] = makeUniformMesh(N, M, type)
bctype = [...
    NdgEdgeType.ZeroGrad, ...
    NdgEdgeType.ZeroGrad, ...
    NdgEdgeType.ZeroGrad, ...
    NdgEdgeType.ZeroGrad];

xlim = [-1, 1];
ylim = [-1, 1];
if (type == NdgCellType.Tri)
    mesh = makeUniformTriMesh(N, xlim, ylim, M, M, bctype);
elseif(type == NdgCellType.Quad)
    mesh = makeUniformQuadMesh(N, xlim, ylim, M, M, bctype);
else
    msgID = [mfile, ':inputCellTypeError'];
    msgtext = 'The input cell type should be NdgCellType.Tri or NdgCellType.Quad.';
    ME = MException(msgID, msgtext);
    throw(ME);
end
end% func

