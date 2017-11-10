classdef ConstAdvSubUniformMesh2d < NdgPhysMat
    %CONSTADVSUBUNIFORMMESH2D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties( Constant )
        Nfield = 1
        Nvar = 1
        varFieldIndex = 1
    end
    
    properties
        %> Number of basis function
        N
        %> Number of elements on each axis
        M
    end
    
    properties( Constant )
        x0 = -0.5;
        y0 = -0.5;
        u0 = 0.5;
        v0 = 0.5;
    end
    
    methods
        function obj = ConstAdvSubUniformMesh2d(N, M)
            mesh = makeSubUniformMesh(N, M);
            obj = obj@NdgPhysMat();
            obj.N = N;
            obj.M = M;
            obj.setNdgPhys( mesh );
        end
    end
    
    methods( Access = protected )
        
        matEvaluateRHS( obj, fphys );
        matEvaluateCellInnerFlux( obj, fphys );
        matEvaluateCellSurfFlux( obj, fphys, fext );
        
        function [ fphys ] = setInitialField( obj )
            fphys = cell( obj.Nmesh, 1 );
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                fphys{m} = getExtFunc(obj, obj.meshUnion(m), 0);
                fphys{m} = mesh.cell.P * fphys{m};
            end
        end% func
        
        function [ option ] = setOption( obj, option )
            outputIntervalNum = 50;
            option('startTime') = 0.0;
            option('finalTime') = 2.0;
            option('LimiterType') = NdgLimiterType.None;
            option('temporalDiscreteType') = NdgIntervalType.Constant;
            option('timeInterval') ...
                = 2/obj.M/sqrt(obj.u0 ^ 2 + obj.v0 ^2)/(2*obj.N + 1);
            option('obcType') = NdgBCType.None;
            option('outputIntervalType') = NdgIntervalType.DeltaTime;
            option('outputTimeInterval') = 2/outputIntervalNum;
            option('outputNetcdfCaseName') = 'ConstAdvUniformMesh2d';
            option('temporalDiscreteType') = NdgTemporalDiscreteType.RK45;
        end
        
        function [E, G] = matEvaluateFlux( obj, mesh, fieldValue )
            E = obj.u0 .* fieldValue;
            G = obj.v0 .* fieldValue;
        end
        
        function [ dflux ] = matEvaluateNumericalFlux( obj, mesh, fphys, fext )
            [ fm ] = fphys(mesh.eidM);
            [ fp ] = fphys(mesh.eidP);
            [ Em, Gm ] = obj.matEvaluateFlux( mesh, fm );
            [ uNorm ] = obj.u0 .* mesh.nx + obj.v0 .* mesh.ny;
            fluxS = ( fm .* ( sign( uNorm ) + 1 ) * 0.5 ...
                + fp .* ( 1 - sign( uNorm )  ) * 0.5 ) .* uNorm;
            dflux = Em .* mesh.nx + Gm .* mesh.ny - fluxS;
        end
        
        %> the exact function
        function [ f_ext ] = getExtFunc(obj, mesh, time)
            xc = obj.x0 + obj.u0.*time;
            yc = obj.y0 + obj.v0.*time;
            
            sigma = 125*1e3/(33*33);
            t = -( (mesh.x-xc).^2+(mesh.y-yc).^2 )*sigma;
            f_ext = exp(t);
        end% func
    end
    
end

function mesh = makeSubUniformMesh(N, M)
bctype = [NdgEdgeType.Clamped, NdgEdgeType.Clamped, ...
    NdgEdgeType.Clamped, NdgEdgeType.Clamped];

% if (type == NdgCellType.Tri)
mesh = makeUniformSubTriMesh(N, [-1, 1], [-1, 1], ...
    M, M, bctype);
% elseif(type == NdgCellType.Quad)
%     mesh = makeUniformQuadMesh(N, [-1, 1], [-1, 1], ...
%         M, M, bctype);
% else
%     msgID = 'ConstAdvSubUniformMesh2d:inputCellTypeError';
%     msgtext = ['The input cell type should be NdgCellType.Tri',...
%         ' or NdgCellType.Tri.'];
%     ME = MException(msgID, msgtext);
%     throw(ME);
% end
end% func

