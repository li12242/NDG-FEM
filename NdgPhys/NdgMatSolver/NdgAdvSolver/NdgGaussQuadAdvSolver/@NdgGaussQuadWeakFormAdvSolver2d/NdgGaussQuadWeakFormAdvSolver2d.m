classdef NdgGaussQuadWeakFormAdvSolver2d < NdgAbstractGaussQuadAdvSolver
    
    methods
        function obj = NdgGaussQuadWeakFormAdvSolver2d( phys )
            obj = obj@NdgAbstractGaussQuadAdvSolver( phys );
            
            for m = 1:phys.Nmesh
                obj.Dr{m} = obj.Dr{m}';
                obj.Ds{m} = obj.Ds{m}';
                obj.Dt{m} = obj.Dt{m}';
            end
        end
        
        function evaluateAdvectionRHS( obj, fphys )
            phys = obj.phys;
            for m = 1:phys.Nmesh % calculate RHS term on each mesh
                mesh = phys.meshUnion(m);
                [ E, G ] = phys.matEvaluateFlux( mesh, fphys{m} );                
                [ fluxS ] = phys.matEvaluateSurfNumFlux( mesh, obj.nx{m}, obj.ny{m}, fphys{m}, phys.fext{m} );
                
                for i = 1:phys.Nvar
                    [ phys.frhs{m}(:,:,i) ] = ...
                        + obj.Dr{m} * ( obj.rxwJ{m}.* (obj.Vq{m} * E(:,:,i)) + obj.rywJ{m}.* ( obj.Vq{m} * G(:,:,i) ) ) ...
                        + obj.Ds{m} * ( obj.sxwJ{m}.* (obj.Vq{m} * E(:,:,i)) + obj.sywJ{m}.* ( obj.Vq{m} * G(:,:,i) ) ) ...
                        - ( obj.LIFT{m} * ( obj.wJs{m} .* ( obj.FVfq{m} * fluxS(:,:,i) ) ));
                    
                    phys.frhs{m}(:,:,i) = sum(bsxfun(@times, obj.invM{m}, phys.frhs{m}(:,:,i)'), 2);
                    
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

