classdef PerturbHump2d < SWEPreBlanaced2d
    
    properties( Constant )
        gra = 9.81
        hmin = 1e-4
    end
    
    methods
        function obj = PerturbHump2d(N, M, cellType)
            [ mesh ] = makeUniformMesh(N, M, cellType);
            obj = obj@SWEPreBlanaced2d();
            obj.initPhysFromOptions( mesh );
            obj.fext = obj.setInitialField( );
        end
    end
    
    methods(Access=protected)
        function fphys = setInitialField( obj )
            fphys = cell( obj.Nmesh, 1 );
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                fphys{m} = zeros( mesh.cell.Np, mesh.K, obj.Nfield );
                
                bot = 0.8*exp( -5*(mesh.x - 0.9).^2 - 50*(mesh.y - 0.5).^2 );
                fphys{m}(:,:,4) = bot;
                fphys{m}(:,:,1) = 1 - bot;
                ind = ( any(mesh.x > 0.05 ) ) & ( any(mesh.x < 0.15 ) );
                fphys{m}(:, ind, 1) = fphys{m}(:, ind, 1) + 0.01;
            end
        end
        
        function [ option ] = setOption( obj, option )
            ftime = 0.6;
            outputIntervalNum = 200;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputCaseName') = mfilename;
            option('temporalDiscreteType') = enumTemporalDiscrete.RK45;
            option('limiterType') = enumLimiter.Vert;
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
            option('CoriolisType') = enumSWECoriolis.None;
            option('WindType') = enumSWEWind.None;
            option('FrictionType') = enumSWEFriction.None;
            option('SWELimiterType') = enumSWELimiter.OnElevation;
        end
    end
end

function [ mesh ] = makeUniformMesh(N, M, type)
bctype = [...
    NdgEdgeType.SlipWall, ...
    NdgEdgeType.SlipWall, ...
    NdgEdgeType.ClampedDepth, ...
    NdgEdgeType.ClampedDepth];

xlim = [0, 2]; 
ylim = [0, 1];
if (type == enumStdCell.Tri)
    mesh = makeUniformTriMesh(N, xlim, ylim, M, ceil(M/2), bctype);
elseif(type == enumStdCell.Quad)
    mesh = makeUniformQuadMesh(N, xlim, ylim, M, ceil(M/2), bctype);
else
    msgID = [mfile, ':inputCellTypeError'];
    msgtext = [ 'The input cell type should be ', ...
        'NdgCellType.Tri or NdgCellType.Quad.' ];
    ME = MException(msgID, msgtext);
    throw(ME);
end
end% func

