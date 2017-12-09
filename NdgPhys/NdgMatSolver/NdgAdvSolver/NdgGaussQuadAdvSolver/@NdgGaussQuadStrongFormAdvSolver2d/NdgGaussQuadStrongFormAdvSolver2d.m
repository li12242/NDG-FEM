classdef NdgGaussQuadStrongFormAdvSolver2d < NdgAbstractGaussQuadAdvSolver
    
    methods
        function obj = NdgGaussQuadStrongFormAdvSolver2d( phys )
            obj = obj@NdgAbstractGaussQuadAdvSolver( phys );
            
            for m = 1:phys.Nmesh
                obj.Vq{m} = obj.Vq{m}';
            end
        end
        
        function evaluateAdvectionRHS( obj, fphys )
            phys = obj.phys;
            for m = 1:phys.Nmesh % calculate RHS term on each mesh
                mesh = phys.meshUnion(m);
                [ E, G ] = phys.matEvaluateFlux( mesh, fphys{m} );
                [ fluxS ] = phys.matEvaluateSurfNumFlux( mesh, obj.nx{m}, obj.ny{m}, fphys{m}, phys.fext{m} );
                [ flux ] = phys.matEvaluateSurfFlux( mesh, obj.nx{m}, obj.ny{m}, fphys{m} );
                
                for i = 1:phys.Nvar
                    [ phys.frhs{m}(:,:,i) ] = ...
                        obj.Vq{m} * (...
                        - obj.rxwJ{m}.*( obj.Dr{m} * E(:,:,i) ) ...
                        - obj.sxwJ{m}.*( obj.Ds{m} * E(:,:,i) ) ...
                        - obj.rywJ{m}.*( obj.Dr{m} * G(:,:,i) ) ...
                        - obj.sywJ{m}.*( obj.Ds{m} * G(:,:,i) ) ) ...
                        + ( obj.LIFT{m} * ( obj.wJs{m} .* ( flux(:,:,i) - fluxS(:,:,i) ) ));
                    
                    %> times the inverse mass matrix to each cell rhs term
%                     for k = 1:mesh.K
%                         phys.frhs{m}(:,k,i) = obj.invM{m}(:,:,k) * phys.frhs{m}(:,k,i);
%                     end
                    phys.frhs{m}(:,:,i) = permute( sum( ...
                        bsxfun(@times, obj.invM{m}, ...
                        permute( permute( phys.frhs{m}(:,:,i), [1,3,2] ), ...
                        [2,1,3] ) ), 2 ), [1,3,2]);
                end
                
            end
        end
    end
    
end

