classdef NdgGaussQuadWeakFormAdvSolver2d < NdgGaussQuadWeakFormSolver & NdgAbstractAdvSolver
    
    methods
        function obj = NdgGaussQuadWeakFormAdvSolver2d( phys )
            obj = obj@NdgGaussQuadWeakFormSolver( phys );
            obj = obj@NdgAbstractAdvSolver( phys );
        end
        
        function evaluateAdvectionRHS( obj, fphys )
            phys = obj.phys;
            for m = 1:phys.Nmesh % calculate RHS term on each mesh
                mesh = phys.meshUnion(m);
                fq = zeros( mesh.cell.Nq, mesh.K, phys.Nfield );
                for i = 1:phys.Nfield
                    fq(:,:,i) = obj.Vq{m} * fphys{m}(:,:,i);
                end
                [ E, G ] = phys.matEvaluateFlux( mesh, fq );
                [ fm, fp ] = phys.matEvaluateSurfaceValue( mesh, fphys{m}, phys.fext{m} );
                fmq = zeros( obj.TNfq{m}, mesh.K, phys.Nfield );
                fpq = zeros( obj.TNfq{m}, mesh.K, phys.Nfield );
                for i = 1:phys.Nvar
                    fmq(:,:,i) = obj.FVfq{m} * fm(:,:,i);
                    fpq(:,:,i) = obj.FVfq{m} * fp(:,:,i);
                end
                [ fluxS ] = phys.matEvaluateSurfNumFlux( mesh, obj.nx{m}, obj.ny{m}, fmq, fpq );
                
                for i = 1:phys.Nvar
                    [ phys.frhs{m}(:,:,i) ] = ...
                        + obj.Dr{m} * ( obj.rxwJ{m}.* (E(:,:,i)) + obj.rywJ{m}.* ( G(:,:,i) ) ) ...
                        + obj.Ds{m} * ( obj.sxwJ{m}.* (E(:,:,i)) + obj.sywJ{m}.* ( G(:,:,i) ) ) ...
                        - ( obj.LIFT{m} * ( obj.wJs{m} .* ( fluxS(:,:,i) ) ));
                    
                    phys.frhs{m}(:,:,i) = permute( sum( ...
                        bsxfun(@times, obj.invM{m}, ...
                        permute( permute( phys.frhs{m}(:,:,i), [1,3,2] ), ...
                        [2,1,3] ) ), 2 ), [1,3,2]);
                end
            end
        end
    end
    
end

