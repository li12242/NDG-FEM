classdef AdvRotationUniformMesh2d < AdvAbstractVarFlow2d
    
    properties
        %> Number of basis function
        N
        %> Number of elements on each axis
        M
    end
    
    properties(Constant)
        %> domain central
        x0 = 0
        %> domain central
        y0 = 0
        %> distance of the initial Gauss mount from the central point
        rd = 0.5
        %> size of the initial Gauss mount
        r0 = 0.25
        %> angle velocity
        w = 5*pi/6;
    end
    
    methods
        function obj = AdvRotationUniformMesh2d(N, M, cellType)
            mesh = makeUniformMesh(N, M, cellType);
            obj = obj@AdvAbstractVarFlow2d( );
            obj.N = N;
            obj.M = M;
            obj.initPhysFromOptions( mesh );
        end
        
        function [ fm, fp ] = matEvaluateEdgeValue( obj, edge, fphys, fext )
            fm = zeros( edge.Nfp, edge.Ne, obj.Nfield );
            fp = zeros( edge.Nfp, edge.Ne, obj.Nfield );
            
            for k = 1:edge.Ne
                for fld = 1:obj.Nfield
                    fm(:, k, fld) = fphys(edge.FToN1(:, k), edge.FToE(1, k), fld);
                    fp(:, k, fld) = fphys(edge.FToN2(:, k), edge.FToE(2, k), fld);
                end
            end
        end
        
        function [ fluxM, fluxP ] = matEvaluateEdgeFlux( obj, edge, fm, fp )
            Em = fm(:,:,1) .* fm(:,:,2);
            Gm = fm(:,:,1) .* fm(:,:,3);
            fluxM = Em .* edge.nx + Gm .* edge.ny;
            
            Em = fp(:,:,1) .* fp(:,:,2);
            Gm = fp(:,:,1) .* fp(:,:,3);
            fluxP = Em .* edge.nx + Gm .* edge.ny;
        end
        
        function [ fluxS ] = matEvaluateEdgeNumFlux( obj, edge, fm, fp )
            [ uNorm ] = fm(:,:,2) .* edge.nx + fm(:,:,3) .* edge.ny;
            sign_um = sign( uNorm );
            fluxS = ( fm(:,:,1).*( sign_um + 1 )*0.5 + fp(:,:,1).*( 1 - sign_um  )*0.5 ).*uNorm;
        end
        
        function [ frhs ] = matEvaluateEdgeRHS( obj, edge, fluxM, fluxP, fluxS )
            frhs = zeros( edge.mesh.cell.Np, edge.mesh.K, obj.Nvar );
            for fld = 1:obj.Nvar
                for k = 1:edge.Ne
                    e1 = edge.FToE(1, k);
                    e2 = edge.FToE(2, k);
                    n1 = edge.FToN1(:, k);
                    n2 = edge.FToN2(:, k);
                    
                    deltaFlux1 = fluxM(:,k,fld) - fluxS(:,k,fld);
                    frhs(n1, e1, fld) = frhs(n1, e1, fld) + ...
                        edge.bcell.M * ( edge.Js(:, k) .* deltaFlux1 );
                    %                     frhs(:, e1, fld) = ...
                    %                         frhs(:, e1, fld) + edge.mesh.cell.invM(:, n1) * ...
                    %                         ( edge.bcell.M * ( edge.Js(:, k) .* deltaFlux1 ) ) ...
                    %                         ./ edge.mesh.J( :, e1 );
                    
                    deltaFlux2 = fluxP(:,k,fld) - fluxS(:,k,fld);
                    frhs(n2, e2, fld) = frhs(n2, e2, fld) - ...
                        edge.bcell.M * ( edge.Js(:, k) .* deltaFlux2 );
                    %                     frhs(:, e2, fld) = ...
                    %                         frhs(:, e2, fld) - edge.mesh.cell.invM(:, n2) * ...
                    %                         ( edge.bcell.M * ( edge.Js(:, k) .* deltaFlux2 ) ) ...
                    %                         ./ edge.mesh.J( :, e2 );
                end
                frhs(:, :, fld) = edge.mesh.cell.invM * frhs(:, :, fld) ./ edge.mesh.J;
            end            
        end
    end
    
    methods( Access = protected )
        
        function [ fphys ] = setInitialField( obj )
            fphys = cell( obj.Nmesh, 1 );
            for m = 1:obj.Nmesh
                fphys{m} = getExtFunc(obj, obj.meshUnion(m), 0);
            end
        end% func
        
        function option = setOption( obj, option )
            outputIntervalNum = 50;
            option('startTime') = 0.0;
            option('finalTime') = 1.2;
            option('obcType') = NdgBCType.None;
            option('outputIntervalType') = NdgIOIntervalType.DeltaTime;
            option('outputTimeInterval') = 2.4/outputIntervalNum;
            option('outputNetcdfCaseName') = mfilename;
            option('temporalDiscreteType') = NdgTemporalDiscreteType.RK45;
            option('timeInterval') = sqrt(2)/obj.M/obj.w/(2*obj.N + 1);
            option('equationType') = NdgDiscreteEquationType.Strong;
            option('integralType') = NdgDiscreteIntegralType.QuadratureFree;
            option('limiterType') = NdgLimiterType.None;
        end
        
        %> the exact function
        function f_ext = getExtFunc(obj, mesh, time)
            f_ext = zeros( mesh.cell.Np, mesh.K, obj.Nfield );
            
            theta0 = -pi;
            theta = theta0 + obj.w*time;
            xt = obj.x0 + obj.rd*cos(theta);
            yt = obj.y0 + obj.rd*sin(theta);
            r2 = sqrt((mesh.x - xt).^2+(mesh.y - yt).^2)./obj.r0;
            ind = ( r2 <= 1.0);
            temp = zeros( mesh.cell.Np, mesh.K );
            temp(ind) = ( 1+cos(r2(ind)*pi) )./2;
            
            f_ext(:,:,1) = temp;
            f_ext(:,:,2) = obj.w.* (- mesh.y);
            f_ext(:,:,3) = obj.w.*( mesh.x );
        end% func
    end% methods
    
end% class

function mesh = makeUniformMesh(N, M, type)
bctype = [...
    NdgEdgeType.Clamped, ...
    NdgEdgeType.Clamped, ...
    NdgEdgeType.Clamped, ...
    NdgEdgeType.Clamped];

if (type == NdgCellType.Tri)
    mesh = makeUniformTriMesh(N, [-1, 1], [-1, 1], M, M, bctype);
elseif(type == NdgCellType.Quad)
    mesh = makeUniformQuadMesh(N, [-1, 1], [-1, 1], M, M, bctype);
else
    msgID = [mfilename, ':inputCellTypeError'];
    msgtext = ['The input cell type should be NdgCellType.Tri',...
        ' or NdgCellType.Tri.'];
    ME = MException(msgID, msgtext);
    throw(ME);
end
end% func

