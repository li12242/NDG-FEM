classdef NdgQuadFreeStrongFormAdvSolver3d < NdgQuadFreeStrongFormSolver & ...
        NdgAbstractAdvSolver
    
    methods
        function obj = NdgQuadFreeStrongFormAdvSolver3d( phys )
            obj = obj@NdgQuadFreeStrongFormSolver( phys );
            obj = obj@NdgAbstractAdvSolver(phys);
        end
        %> Call the flux subroutine from the NdgPhys object.
        function evaluateAdvectionRHS( obj, fphys )
            phys = obj.phys;
            
            % evaluate inner edge
            for m = 1:phys.Nmesh
                edge = phys.meshUnion(m).InnerEdge;
                [ fm, fp ] = edge.matEvaluateSurfValue( fphys );
                [ fluxM ] = phys.matEvaluateSurfFlux( edge, edge.nx, edge.ny, edge.nz, fm );
                [ fluxP ] = phys.matEvaluateSurfFlux( edge, edge.nx, edge.ny, edge.nz, fp );
                [ fluxS ] = phys.matEvaluateSurfNumFlux( edge, edge.nx, edge.ny, edge.nz, fm, fp );
                [ phys.frhs{m} ] = edge.matEvaluateStrongFromEdgeRHS( fluxM, fluxP, fluxS );

                edge = phys.meshUnion(m).BottomEdge;
                [ fm, fp ] = edge.matEvaluateSurfValue( fphys );
                %[ fm, fp ] = phys.matImposeBoundaryCondition( edge, edge.nx, edge.ny, fm, fp, phys.fext );
                [ fluxM ] = phys.matEvaluateSurfFlux( edge, edge.nx, edge.ny, edge.nz, fm );
                [ fluxP ] = phys.matEvaluateSurfFlux( edge, edge.nx, edge.ny, edge.nz, fp );
                [ fluxS ] = phys.matEvaluateSurfNumFlux( edge, edge.nx, edge.ny, edge.nz, fm, fp );
                [ phys.frhs{m} ] = phys.frhs{m} + edge.matEvaluateStrongFormEdgeRHS( fluxM, fluxP, fluxS );
            end
            
            for m = 1:phys.Nmesh % calculate RHS term on each mesh
                mesh = phys.meshUnion(m);
                [ E, G, H ] = phys.matEvaluateFlux( mesh, fphys{m} );
                
                for i = 1:phys.Nvar
                    phys.frhs{m}(:,:,i) = ...
                        phys.frhs{m}(:,:,i) + ...
                        - obj.rx{m}.*( obj.Dr{m} * E(:,:,i) ) ...
                        - obj.sx{m}.*( obj.Ds{m} * E(:,:,i) ) ...
                        - obj.tx{m}.*( obj.Dt{m} * E(:,:,i) ) ...
                        - obj.ry{m}.*( obj.Dr{m} * G(:,:,i) ) ...
                        - obj.sy{m}.*( obj.Ds{m} * G(:,:,i) ) ...
                        - obj.ty{m}.*( obj.Dt{m} * G(:,:,i) ) ...
                        - obj.rz{m}.*( obj.Dr{m} * H(:,:,i) ) ...
                        - obj.sz{m}.*( obj.Ds{m} * H(:,:,i) ) ...
                        - obj.tz{m}.*( obj.Dt{m} * H(:,:,i) );
                end
            end
        end
    end
    
end