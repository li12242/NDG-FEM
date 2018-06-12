
classdef ConicalLandRunup2d < SWEWDPreBlanaced2d
    
    properties(Constant)
        %> wet/dry depth threshold
        hmin = 1e-4
        %> gravity acceleration
        gra = 9.8
        %> initial water depth
        h0 = 0.32
        %> island central
        x0 = 12.5
        y0 = 15
        %> island length
        L = 15
        alpha = 0.1
    end
    
    methods
        function obj = ConicalLandRunup2d(N, M, cellType)
            [ mesh ] = makeUniformMesh(N, M, cellType);
            % mesh = readGmshFile( N );
            obj = obj@SWEWDPreBlanaced2d();
            obj.initPhysFromOptions( mesh );
        end
        
        drawMaxRunup( obj )
        drawResult( obj )
        drawGaugeResult( obj )
    end
    
    methods(Access=protected)
        function matUpdateExternalField( obj, time, fphys )
            beta = ( obj.h0/obj.L )^2;
            tx = sqrt( 3 * obj.alpha / 4 / beta * ( 1 + obj.alpha ) ) ...
                * sqrt( obj.gra * obj.h0 ) / obj.L * ( time );
%             alpha = obj.alpha;
%             L = obj.L;
%             h0 = obj.h0;
%             xi = sqrt( 3 * alpha * ( 1+alpha ) * L^2 / 4 / h0^2 );
%             tx = xi * sqrt( obj.gra * h0 / L ) * ( time );
            
            sech2 = ( 2 / ( exp(tx) + exp(-tx) ) )^2;
            theta = 1.9;
            h = obj.h0 + obj.alpha * obj.h0 * sech2 * theta;
            for m = 1:obj.Nmesh
                obj.fext{m}(:,:,1) = h;
%                 obj.fext{m}(:,:,1) = obj.h0 - obj.fext{m}(:,:,4);
            end
        end
        
        function fphys = setInitialField( obj )
            fphys = cell( obj.Nmesh, 1 );
            xt = obj.x0;
            yt = obj.y0;
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                fphys{m} = zeros( mesh.cell.Np, mesh.K, obj.Nfield );
                
                r = sqrt( (mesh.x - xt).^2 + (mesh.y - yt).^2 );
                bot = min( 0.625, 0.9 - r/4 );
                bot( bot <= 0 ) = 0.0;
                h = obj.h0 - bot;
                h( h < 0 ) = 0.0;
                
                fphys{m}(:,:,1) = h;
                fphys{m}(:,:,4) = bot;
            end
        end
        
        function [ option ] = setOption( obj, option )
            ftime = 20;
            outputIntervalNum = 100;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('outputIntervalType') = enumOutputInterval.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputCaseName') = mfilename;
            option('temporalDiscreteType') = enumTemporalDiscrete.RK45;
            option('limiterType') = enumLimiter.Vert;
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
            option('SWELimiterType') = enumSWELimiter.OnElevation;
        end
    end
    
end

function [ mesh ] = readGmshFile( N )
[ path, name, ext ] = fileparts( mfilename('fullpath') );
obj.gmshFile = [ path, '/mesh/conicalLand.msh' ];
mesh = makeGmshFileUMeshUnion2d( N, obj.gmshFile );
end

function [ mesh ] = makeUniformMesh(N, M, type)
bctype = [...
    enumBoundaryCondition.ZeroGrad, ...
    enumBoundaryCondition.ZeroGrad, ...
    enumBoundaryCondition.ClampedDepth, ...
    enumBoundaryCondition.ZeroGrad];

if (type == enumStdCell.Tri)
    mesh = makeUniformTriMesh(N, [0, 25], [0, 30], M, M, bctype);
elseif(type == enumStdCell.Quad)
    mesh = makeUniformQuadMesh(N, [0, 25], [0, 30], M, M, bctype);
else
    msgID = [mfile, ':inputCellTypeError'];
    msgtext = 'The input cell type should be NdgCellType.Tri or NdgCellType.Quad.';
    ME = MException(msgID, msgtext);
    throw(ME);
end
end% func

