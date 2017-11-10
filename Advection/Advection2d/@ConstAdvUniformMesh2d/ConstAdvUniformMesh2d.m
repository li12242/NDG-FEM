%> The advection problem with the contant flow filed.
%> 
classdef ConstAdvUniformMesh2d < NdgPhysMat
    
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
        function obj = ConstAdvUniformMesh2d(N, M, cellType)
            
            mesh = makeUniformMesh(N, M, cellType);
            obj = obj@NdgPhysMat();
            obj.N = N;
            obj.M = M;
            obj.setNdgPhys( mesh );            
        end
        
        function drawResult(obj, n)
            filename = [obj.getOption('outputNetcdfCaseName'), '.1-1.nc'];
            time = ncread(filename, 'time');
            Ntime = numel(time);

            for i = 1:n:Ntime
                field = obj.accessOutputResultAtStepNum( i );
                obj.meshUnion.draw( field{1}(:,:,1) ); 
                drawnow;
            end
        end% func

    end
    
    methods( Access = protected )
        function [ fphys ] = setInitialField( obj )
            fphys = cell( obj.Nmesh, 1 );
            for m = 1:obj.Nmesh
                fphys{m} = getExtFunc(obj, obj.meshUnion(m), 0);
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
        
        function matEvaluateRHS( obj, fphys )
            obj.matEvaluateRHS2d( fphys );
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

function mesh = makeUniformMesh(N, M, type)
bctype = [NdgEdgeType.Clamped, NdgEdgeType.Clamped, ...
    NdgEdgeType.Clamped, NdgEdgeType.Clamped];

if (type == NdgCellType.Tri)
    mesh = makeUniformTriMesh(N, [-1, 1], [-1, 1], ...
        M, M, bctype);
elseif(type == NdgCellType.Quad)
    mesh = makeUniformQuadMesh(N, [-1, 1], [-1, 1], ...
        M, M, bctype);
else
    msgID = 'AdvUniformUnion:inputCellTypeError';
    msgtext = ['The input cell type should be NdgCellType.Tri',...
        ' or NdgCellType.Tri.'];
    ME = MException(msgID, msgtext);
    throw(ME);
end
end% func

