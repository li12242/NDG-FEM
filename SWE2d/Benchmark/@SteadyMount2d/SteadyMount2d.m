classdef SteadyMount2d < SWEPreBlanaced2d %& SDBAbstractTest & CSBAbstractTest
    
    properties( Constant )
        %> wet/dry depth threshold
        hmin = 1e-2
        %> gravity acceleration
        gra = 9.8
        
        x0 = 0.5
        y0 = 0.5
    end
    
    methods
        function obj = SteadyMount2d( N, M, cellType )
            [ mesh ] = makeUniformMesh(N, M, cellType);
            obj = obj@SWEPreBlanaced2d();
            obj.initPhysFromOptions( mesh );
            obj.fphys = obj.matEvaluatePostFunc( obj.fphys );
            obj.fext = obj.setBoundaryCondition( );
        end
        
        
    end
    methods(Access=protected)
        function fext = setBoundaryCondition( obj )
            fext = cell( obj.Nmesh, 1 );
            for m = obj.Nmesh
                mesh = obj.meshUnion(m);
                fext{m} = zeros( mesh.cell.Np, mesh.K, obj.Nfield );
                fext{m}(:, :, 1) = 0.2;
                
            end
        end
        
        function fphys = setInitialField( obj )
            fphys = getExactFunction(obj, 0);
        end
        
        function [ option ] = setOption( obj, option )
            ftime = 10;
            outputIntervalNum = 50;
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
        
%         function [ fphys ] = matEvaluatePostFunc( obj, fphys )
%             for m = 1:obj.Nmesh
%                 hc = obj.meshUnion(m).GetMeshAverageValue( fphys{m}(:,:,1) );
%                 qxc = obj.meshUnion(m).GetMeshAverageValue( fphys{m}(:,:,2) );
%                 qyc = obj.meshUnion(m).GetMeshAverageValue( fphys{m}(:,:,3) );
%                 fphys{m}(:,:,1:3) = mxEvaluatePostFunc2d( obj.hmin, fphys{m}, hc, qxc, qyc );
%                 fphys{m}(:,:,2) = fphys{m}(:,:,2) .* 0.9;
%                 fphys{m}(:,:,3) = fphys{m}(:,:,3) .* 0.9;
%             end
%             obj.matUpdateWetDryState( fphys );
%         end
        
        function fphys = getExactFunction( obj, time )
            
            fphys = cell( obj.Nmesh, 1 );
%             limiter = NdgVertLimiter2d( obj.meshUnion );
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                bot = 0.25 - 2.5*(mesh.x - obj.x0).^2 - 2.5*(mesh.y - obj.y0).^2;
                bmin = min(bot);
                bmean = mesh.GetMeshAverageValue( bot );
                bmean( bmean <= 0 ) = 0;
                theta = bsxfun(@min, 1, bmean./(bmean - bmin) );
                bot = bsxfun(@plus, bsxfun(@times, bsxfun(@plus, bot, -bmean ), theta ), bmean);
% %                 sigma = 25*1e3/(33*33);
% %                 temp = sigma*( (mesh.x - obj.x0).^2 + (mesh.y - obj.y0).^2 );
% %                 bot = exp( -temp );
%                 fphys{m}(:,:,4) = bot;
%             end
%             fphys = limiter.matLimit( fphys, 4 );
%             for m = 1:obj.Nmesh
%                 bot = fphys{m}(:,:,4);
                eta = 0.2 * ones( mesh.cell.Np, mesh.K );
                h = max( eta - bot, 0 );
                ind = ( max(h) > 0 );
                bot(:, ind) = eta(:, ind) - h(:, ind);
                fphys{m}(:,:,1) = h;
                fphys{m}(:,:,4) = bot;
            end
        end
    end
    
end

function [ mesh ] = makeUniformMesh(N, M, type)
bctype = [...
    NdgEdgeType.SlipWall, ...
    NdgEdgeType.SlipWall, ...
    NdgEdgeType.SlipWall, ...
    NdgEdgeType.SlipWall];

if (type == NdgCellType.Tri)
    mesh = makeUniformTriMesh(N, [0, 1], [0, 1], M, M, bctype);
elseif(type == NdgCellType.Quad)
    mesh = makeUniformQuadMesh(N, [0, 1], [0, 1], M, M, bctype);
else
    msgID = [mfile, ':inputCellTypeError'];
    msgtext = 'The input cell type should be NdgCellType.Tri or NdgCellType.Quad.';
    ME = MException(msgID, msgtext);
    throw(ME);
end
end% func
