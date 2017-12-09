classdef NdgQuadFreeStrongAdvSolver1d < NdgAbstractAdvectionSolver
        
    methods
        %> Call the flux subroutine from the NdgPhys object.
        function evaluateAdvectionRHS( obj, physClass, fphys )
            for m = 1:physClass.Nmesh % calculate RHS term on each mesh
                mesh = physClass.meshUnion(m);
                [ E ] = physClass.matEvaluateFlux( mesh, fphys{m} );
                [ dFlux ] = physClass.matEvaluateNumericalFlux( mesh, fphys{m}, physClass.fext{m} );
                
                for i = 1:physClass.Nvar
                    [ physClass.frhs{m}(:,:,i) ] = ...
                        - mesh.rx.*( mesh.cell.Dr * E(:,:,i) ) ...
                        + ( mesh.cell.LIFT * ( mesh.Js .* dFlux(:,:,i) ))./ mesh.J;
                end
                
            end
        end
    end
    
end

