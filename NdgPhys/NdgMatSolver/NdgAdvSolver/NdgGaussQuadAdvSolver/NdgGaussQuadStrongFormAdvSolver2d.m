classdef NdgGaussQuadStrongFormAdvSolver2d < NdgGaussQuadStrongFormSolver & NdgAbstractAdvSolver
    
    methods
        function obj = NdgGaussQuadStrongFormAdvSolver2d( phys )
            obj = obj@NdgGaussQuadStrongFormSolver( phys );
            obj = obj@NdgAbstractAdvSolver( phys );
        end
        
        function evaluateAdvectionRHS( obj, fphys )
            phys = obj.phys;
            for m = 1:phys.Nmesh % calculate RHS term on each mesh
                mesh = phys.meshUnion(m);
                [ E, G ] = phys.matEvaluateFlux( mesh, fphys{m} );
                [ fm, fp ] = phys.matEvaluateSurfaceValue( mesh, fphys{m}, phys.fext{m} );
                fmq = zeros( obj.TNfq{m}, mesh.K, phys.Nfield );
                fpq = zeros( obj.TNfq{m}, mesh.K, phys.Nfield );
                for i = 1:phys.Nvar
                    fmq(:,:,i) = obj.FVfq{m} * fm(:,:,i);
                    fpq(:,:,i) = obj.FVfq{m} * fp(:,:,i);
                end
                [ fluxS ] = phys.matEvaluateSurfNumFlux( mesh, obj.nx{m}, obj.ny{m}, fmq, fpq );
                [ flux ] = phys.matEvaluateSurfFlux( mesh, obj.nx{m}, obj.ny{m}, fmq );
                 
                for i = 1:phys.Nvar
                    [ phys.frhs{m}(:,:,i) ] = ...
                        obj.Vq{m} * (...
                        - obj.rxwJ{m}.*( obj.Dr{m} * E(:,:,i) ) ...
                        - obj.sxwJ{m}.*( obj.Ds{m} * E(:,:,i) ) ...
                        - obj.rywJ{m}.*( obj.Dr{m} * G(:,:,i) ) ...
                        - obj.sywJ{m}.*( obj.Ds{m} * G(:,:,i) ) ) ...
                        + ( obj.LIFT{m} * ( obj.wJs{m} .* ( flux(:,:,i) - fluxS(:,:,i) ) ));
                    
                    phys.frhs{m}(:,:,i) = permute( sum( ...
                        bsxfun(@times, obj.invM{m}, ...
                        permute( permute( phys.frhs{m}(:,:,i), [1,3,2] ), ...
                        [2,1,3] ) ), 2 ), [1,3,2]);
                end
                
            end
        end
    end
    
end
