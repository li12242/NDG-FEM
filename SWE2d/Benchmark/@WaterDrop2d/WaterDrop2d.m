classdef WaterDrop2d < SWEAbstractDBN62d
    %WATERDROP2D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
        gra = 9.81;
        hmin = 1e-4;
    end
    
    methods
        function obj = WaterDrop2d(N, M, cellType)
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
                fphys{m}(:,:,1) = 2.4*(1 + ...
                    exp( -((mesh.x - 10).^2 + (mesh.y - 10).^2)/4 ));
                
            end
        end
        
        function [ option ] = setOption( obj, option )
            ftime = 1.2;
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
        end
    end
    
end

function [ mesh ] = makeUniformMesh(N, M, type)
bctype = [...
    NdgEdgeType.ZeroGrad, ...
    NdgEdgeType.ZeroGrad, ...
    NdgEdgeType.ZeroGrad, ...
    NdgEdgeType.ZeroGrad];

xlim = [0, 20];
ylim = [0, 20];
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

