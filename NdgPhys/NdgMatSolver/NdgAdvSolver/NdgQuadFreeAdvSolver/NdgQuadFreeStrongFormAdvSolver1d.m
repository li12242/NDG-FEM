classdef NdgQuadFreeStrongFormAdvSolver1d < NdgAbstractAdvSolver ...
        & NdgQuadFreeStrongFormSolver
        
    methods
        function obj = NdgQuadFreeStrongFormAdvSolver1d( phys )
            obj = obj@NdgQuadFreeStrongFormSolver( phys );
            obj = obj@NdgAbstractAdvSolver(phys);
        end
        
        %> Call the flux subroutine from the NdgPhys object.
        function evaluateAdvectionRHS( obj, fphys )
            phys = obj.phys;
            for m = 1:phys.Nmesh % calculate RHS term on each mesh
                mesh = phys.meshUnion(m);
                [ E ] = phys.matEvaluateFlux( mesh, fphys{m} );
                [ fm, fp ] = phys.matEvaluateSurfaceValue( mesh, fphys{m}, phys.fext{m} );
                [ fluxS ] = phys.matEvaluateSurfNumFlux( mesh, obj.nx{m}, fm, fp );
                [ flux ] = phys.matEvaluateSurfFlux( mesh, obj.nx{m}, fm );
                
                for i = 1:phys.Nvar
                    [ phys.frhs{m}(:,:,i) ] = ...
                        - mesh.rx.*( mesh.cell.Dr * E(:,:,i) ) ...
                        + ( obj.LIFT{m} * ( obj.Js{m} .* ( flux(:,:,i) - fluxS(:,:,i) ) ))./ obj.J{m};
                end
                
            end
        end
    end
    
end

