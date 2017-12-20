classdef NdgQuadFreeWeakFormAdvSolver2d < NdgQuadFreeWeakFormSlver & ...
        NdgAbstractAdvSolver
    
    methods
        function obj = NdgQuadFreeWeakFormAdvSolver2d( phys )
            obj = obj@NdgQuadFreeWeakFormSlver( phys );
            obj = obj@NdgAbstractAdvSolver( phys );
        end
        
        %> Call the flux subroutine from the NdgPhys object.
        function evaluateAdvectionRHS( obj, fphys )
            phys = obj.phys;
            for m = 1:phys.Nmesh % calculate RHS term on each mesh
                mesh = phys.meshUnion(m);
                [ E, G ] = phys.matEvaluateFlux( mesh, fphys{m} );
                [ fm, fp ] = phys.matEvaluateSurfaceValue( mesh, fphys{m}, phys.fext{m} );
                [ fluxS ] = phys.matEvaluateSurfNumFlux( mesh, obj.nx{m}, obj.ny{m}, fm, fp );
                
                for i = 1:phys.Nvar
                    [ phys.frhs{m}(:,:,i) ] = ...
                        + obj.Dr{m} * ( obj.rx{m}.* (obj.M{m} * E(:,:,i)) + obj.ry{m}.* ( obj.M{m} * G(:,:,i) ) ) ...
                        + obj.Ds{m} * ( obj.sx{m}.* (obj.M{m} * E(:,:,i)) + obj.sy{m}.* ( obj.M{m} * G(:,:,i) ) ) ...
                        - ( obj.LIFT{m} * ( obj.Js{m} .* ( fluxS(:,:,i) ) ))./ obj.J{m};
                end
            end
        end
    end
    
end

