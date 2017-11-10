classdef AdvMixMesh2d < NdgPhysMat

    properties(Constant)
        %> Number of physical field
        Nfield = 3
        %> Number of variable field
        Nvar = 1
        %> field index of variable field
        varFieldIndex = 1
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
        %> gmsh file
        gmshFile = 'Advection/Advection2d/@AdvMixMesh2d/mesh/MixMesh.msh';
    end
    
    properties( SetAccess = protected )
        N
    end
    
    methods
        function obj = AdvMixMesh2d( N )
            obj = obj@NdgPhysMat();
            mesh = makeGmshFileUMeshUnion2d( N, obj.gmshFile );
            obj.N = N;
            obj.setNdgPhys( mesh );
        end
        
        function err = getNormErr2( obj )
            ftime = obj.getOption('finalTime');
            for m = 1:obj.Nmesh
                obj.fext{m} = getExtFunc(obj, obj.meshUnion(m), ftime);
            end
            err = obj.evaluateNormErr2();
        end% func
                
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
        
        function option = setOption( obj, option )
            outputIntervalNum = 50;
            option('startTime') = 0.0;
            option('finalTime') = 2;
            option('LimiterType') = NdgLimiterType.None;
            option('temporalDiscreteType') = NdgIntervalType.Constant;
            option('obcType') = NdgBCType.None;
            option('outputIntervalType') = NdgIntervalType.DeltaTime;
            option('outputTimeInterval') = 2.4/outputIntervalNum;
            option('outputNetcdfCaseName') = 'AdvUniformMesh2d';
            option('temporalDiscreteType') = NdgTemporalDiscreteType.RK45;
            
            dx = 2;
            for m = 1:obj.Nmesh
                dx = min( dx, min( obj.meshUnion(m).charLength ) );
            end
            option('timeInterval') = dx/sqrt(2)/obj.w/(2*obj.N + 1)/2;
        end
        
        function matEvaluateRHS( obj, fphys )
            obj.matEvaluateRHS2d( fphys );
        end
        
        function [E, G] = matEvaluateFlux( obj, mesh, fphys )
            E = fphys(:,:,2) .* fphys(:,:,1);
            G = fphys(:,:,3) .* fphys(:,:,1);
        end
                
        function [dflux] = matEvaluateNumericalFlux( obj, mesh, fphys, fext )
            Ntp = mesh.cell.Np * mesh.K;
            [ fm ] = fphys(mesh.eidM); 
            [ fp ] = fphys(mesh.eidP);
            [ fpext ] = fext(mesh.eidP);
            ind = ( mesh.eidtype == int8(NdgEdgeType.GaussEdge) );
            fp(ind) = fpext(ind);
            [ um ] = fphys(mesh.eidM + Ntp); 
            [ vm ] = fphys(mesh.eidM + 2*Ntp); 
            Em = fm .* um;
            Gm = fm .* vm;
            [ uNorm ] = um .* mesh.nx + vm .* mesh.ny;
            fluxS = ( fm .* ( sign( uNorm ) + 1 ) * 0.5 ...
                + fp .* ( 1 - sign( uNorm )  ) * 0.5 ) .* uNorm;
            dflux = Em .* mesh.nx + Gm .* mesh.ny - fluxS;
        end
        
        %> the exact function
        function f_ext = getExtFunc(obj, mesh, time)
            xc = -0.5 + 0.5*time;
            yc = -0.5 + 0.5*time;
            r2 = sqrt((mesh.x - xc).^2+(mesh.y - yc).^2)./obj.r0;
            ind = ( r2 <= 1.0);
            temp = zeros( mesh.cell.Np, mesh.K );
            temp(ind) = ( 1+cos(r2(ind)*pi) )./2;
            
            f_ext(:,:,1) = temp;
            f_ext(:,:,2) = 0.5;
            f_ext(:,:,3) = 0.5;
        end
%         function f_ext = getExtFunc(obj, mesh, time)
%             f_ext = zeros( mesh.cell.Np, mesh.K, obj.Nfield );
%             
%             theta0 = -pi/2;
%             theta = theta0 + obj.w*time;
%             xt = obj.x0 + obj.rd*cos(theta);
%             yt = obj.y0 + obj.rd*sin(theta);
%             r2 = sqrt((mesh.x - xt).^2+(mesh.y - yt).^2)./obj.r0;
%             ind = ( r2 <= 1.0);
%             temp = zeros( mesh.cell.Np, mesh.K );
%             temp(ind) = ( 1+cos(r2(ind)*pi) )./2;
%             
%             f_ext(:,:,1) = temp;
%             f_ext(:,:,2) = obj.w.* (- mesh.y); 
%             f_ext(:,:,3) = obj.w.*( mesh.x );
%         end% func
    end% methods
    
end

