classdef ClosedChannelSponge2d < ClosedChannel2d
    
    properties
        spgLength = 3e3; %> sponge region size
        distance %> distance to boundary nodes
        sigma %> sponge strength
        maxSigma %> maximum sponge strength
    end
    
    methods( Access = protected, Static )
        function obtype = setOpenBoundaryCondition(  )
            obtype = [ NdgEdgeType.ZeroGrad, NdgEdgeType.SlipWall ];
        end
    end
    
    methods( Access = public )
        function obj = ClosedChannelSponge2d( N, M )
            obj = obj@ClosedChannel2d( N, M );
            
            bp = - obj.ChLength + obj.spgLength;
            ind = obj.meshUnion.xc < bp; % left part is sponge region
            obj.meshUnion.EToR(ind) = NdgRegionType.Sponge;
            
            Nb = 10;
            xb = bp * ones( Nb, 1 ); 
            yb = linspace( 0, obj.ChWidth, Nb )';
            obj.evaluateSpongeDistance( xb, yb );
            
            dt = matUpdateTimeInterval( obj, obj.fphys );
            obj.evaluateSpongeStrength( obj.spgLength, 0.9/dt );
        end
        
    end
    
    methods (Access = protected)
        function [ option ] = setOption( obj, option )
            ftime = 5000;
            outputIntervalNum = 2000;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('obcType') = NdgBCType.None;
            option('outputIntervalType') = NdgIOIntervalType.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputNetcdfCaseName') = ...
                [mfilename, '.', num2str(obj.meshUnion.cell.N)];
            option('temporalDiscreteType') = NdgTemporalDiscreteType.RK22;
            option('limiterType') = NdgLimiterType.None;
            option('equationType') = NdgDiscreteEquationType.Strong;
            option('integralType') = NdgDiscreteIntegralType.GaussQuadrature;
            option('CoriolisType') = SWECoriolisType.None;
            option('WindType') = SWEWindType.None;
            option('FrictionType') = SWEFrictionType.None;
        end
        
        function matEvaluateTopographySourceTerm( obj, fphys )
            matEvaluateTopographySourceTerm@ClosedChannel2d( obj, fphys );
            
            for m = 1:obj.Nmesh
                obj.frhs{m}(:,:,1) = obj.frhs{m}(:,:,1)...
                    - obj.sigma.*( fphys{m}(:,:,1) - obj.fext{m}(:,:,1) );
                obj.frhs{m}(:,:,2) = obj.frhs{m}(:,:,2)...
                    - obj.sigma.*( fphys{m}(:,:,2) - obj.fext{m}(:,:,2) );
                obj.frhs{m}(:,:,3) = obj.frhs{m}(:,:,3)...
                    - obj.sigma.*( fphys{m}(:,:,3) - obj.fext{m}(:,:,3) );
            end
        end
        
        %> \brief calculate distance from the boundary
        function evaluateSpongeDistance(obj, xb, yb)
            
            mesh = obj.meshUnion;
            obj.distance = zeros( mesh.cell.Np, mesh.K);
            for k = 1:mesh.K
                if (mesh.EToR(k) ~= NdgRegionType.Sponge)
                    continue;
                end
                
                for n = 1:mesh.cell.Np
                    xi = mesh.x(n, k);
                    yi = mesh.y(n, k);
                    obj.distance(n, k) = min( sqrt( (xi - xb).^2 + (yi - yb).^2 ) );
                end
            end
        end% func
        
        %> \brief
        %> evaluate sponge strength, care may be needed numerically to guarantee that
        %> maxSigma * dt <= 0.9;
        function evaluateSpongeStrength(obj, spongeLength, maxSigma)
            mesh = obj.meshUnion;
            obj.sigma = zeros(mesh.cell.Np, mesh.K);
            p = 3;
            for k = 1:mesh.K
                if mesh.EToR(k) ~= NdgRegionType.Sponge
                    continue;
                end
                obj.sigma(:, k) = maxSigma*abs( obj.distance(:, k)/spongeLength ).^p;
            end
        end% func
    end
    
end

