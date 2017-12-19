%> @brief Dam break problem over a plan.
%
%> This problem consider the positivity-perserving methods of the dam break
%> problem over the non-flat bottom. The exact wet/dry front and its
%> velocity is given in [1] with
%> 
%> [1] Xing Y, Zhang X, Shu C-W. Positivity-preserving high order
%> well-balanced discontinuous Galerkin methods for the shallow water
%> equations. Advances in Water Resources 2010;33:1476–93.
%
%
% ======================================================================
%> This class is part of the NDGOM software. 
%> @author li12242, Tianjin University, li12242@tju.edu.cn
% ======================================================================
classdef DamBreakPlanUmiformMesh2d < SWEAbstractDBN62d
    
    properties(Constant)
        %> wet/dry depth threshold
        hmin = 1e-4
        %> gravity acceleration
        gra = 9.8
        %> Dam position
        damPosition = 0
        %> topography angle
        alpha = pi/60
    end
    
    methods
        function obj = DamBreakPlanUmiformMesh2d(N, M, cellType)
            [ mesh ] = makeUniformMesh(N, M, cellType);
            obj = obj@SWEAbstractDBN62d();
            obj.initPhysFromOptions( mesh );
        end
        
        
    end
    methods(Access=protected)
        function fphys = setInitialField( obj )
            fphys = getExactFunction(obj, 0);
        end
        
        function [ option ] = setOption( obj, option )
            ftime = 2;
            outputIntervalNum = 50;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('temporalDiscreteType') = NdgTemporalIntervalType.DeltaTime;
            option('obcType') = NdgBCType.None;
            option('outputIntervalType') = NdgIOIntervalType.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputNetcdfCaseName') = mfilename;
            option('temporalDiscreteType') = NdgTemporalDiscreteType.RK45;
            option('limiterType') = NdgLimiterType.Vert;
            option('equationType') = NdgDiscreteEquationType.Strong;
            option('integralType') = NdgDiscreteIntegralType.QuadratureFree;
            option('CoriolisType')=CoriolisType.None;
            option('WindType')=WindType.None;
            option('FrictionType')=FrictionType.None;
%             option('WellBlancedType') = true;
        end
        
        function fphys = getExactFunction( obj, time )
            
            fphys = cell( obj.Nmesh, 1 );
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                fphys{m} = zeros( mesh.cell.Np, mesh.K, obj.Nfield );
                bot = mesh.x * tan( obj.alpha );
                fphys{m}(:,:,4) = bot;
                ind = ( mesh.x < 0 );
                temp = zeros( mesh.cell.Np, mesh.K );
                temp( ind ) = 1 - bot( ind );
                fphys{m}(:,:,1) = temp;
            end
        end
    end
    
end

function [ mesh ] = makeUniformMesh(N, M, type)
bctype = [...
    NdgEdgeType.ZeroGrad, ...
    NdgEdgeType.ZeroGrad, ...
    NdgEdgeType.SlipWall, ...
    NdgEdgeType.ZeroGrad];

if (type == NdgCellType.Tri)
    mesh = makeUniformTriMesh(N, [-15, 15], [-15, 15], M, ceil(M/30), bctype);
elseif(type == NdgCellType.Quad)
    mesh = makeUniformQuadMesh(N, [-15, 15], [-15, 15], M, ceil(M/30), bctype);
else
    msgID = [mfile, ':inputCellTypeError'];
    msgtext = 'The input cell type should be NdgCellType.Tri or NdgCellType.Quad.';
    ME = MException(msgID, msgtext);
    throw(ME);
end
end% func

