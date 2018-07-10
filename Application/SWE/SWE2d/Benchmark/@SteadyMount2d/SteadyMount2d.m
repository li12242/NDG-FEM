classdef SteadyMount2d < SWEConventional2d
    
    properties( Constant )
        %> wet/dry depth threshold
        hmin = 1e-2
        %> gravity acceleration
        gra = 9.8
        
        x0 = 0.5
        y0 = 0.5
    end
    
    properties
        N
        M
    end
    
    methods
        function obj = SteadyMount2d( N, M, cellType )
            [ mesh ] = makeUniformMesh(N, M, cellType);
            obj = obj@SWEConventional2d();
            obj.N = N;
            obj.M = M;
            obj.initPhysFromOptions( mesh );
            obj.fphys = obj.matEvaluatePostFunc( obj.fphys );
        end
        
        DrawFluxError( obj );
        CoutourFlux( obj );
        DrawFluxSectionError( obj );
    end
    
    methods(Access=protected)
        function [ fphys ] = matEvaluatePostFunc( obj, fphys )
            fphys{1}(:, :, 2) = fphys{1}(:, :, 2) .* 9e-1;
            fphys{1}(:, :, 3) = fphys{1}(:, :, 3) .* 9e-1;
        end
        
        function fphys = setInitialField( obj )
            fphys = getExactFunction(obj, 0);
        end
        
        function [ option ] = setOption( obj, option )
            ftime = 10;
            outputIntervalNum = 200;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('outputIntervalType') = enumOutputInterval.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputCaseName') = [mfilename, '_', num2str(obj.N)];
            option('temporalDiscreteType') = enumTemporalDiscrete.RK45;
            option('limiterType') = enumLimiter.None;
            option('SWELimiterType') = enumSWELimiter.OnElevation;
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
            option('CoriolisType')= enumSWECoriolis.None;
            option('WindType') = enumSWEWind.None;
            option('FrictionType') = enumSWEFriction.None;
            option('NumFluxType') = enumSWENumFlux.HLL;
        end
        
        function fphys = getExactFunction( obj, time )
            fphys = cell( obj.Nmesh, 1 );
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                bot = 0.25 * exp( 40 * (-(mesh.x - .5).^2 - (mesh.y - 0.5).^2 ));
                eta = ones( mesh.cell.Np, mesh.K );
                
                fphys{m}(:,:,1) = eta - bot;
                fphys{m}(:,:,4) = bot;
            end
        end
    end
    
end

function [ mesh ] = makeUniformMesh(N, M, type)
bctype = [...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall];

if (type == enumStdCell.Tri)
    mesh = makeUniformTriMesh(N, [0, 1], [0, 1], M, M, bctype);
elseif(type == enumStdCell.Quad)
    mesh = makeUniformQuadMesh(N, [0, 1], [0, 1], M, M, bctype);
else
    msgID = [mfile, ':inputCellTypeError'];
    msgtext = 'The input cell type should be NdgCellType.Tri or NdgCellType.Quad.';
    ME = MException(msgID, msgtext);
    throw(ME);
end
end% func
