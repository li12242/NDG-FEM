%> @brief Evaluating the RHS term for the 2d problem
%> @details The function calculate the RHS term on each mesh
function matEvaluateRHS( obj, fphys )

for m = 1:obj.Nmesh % calculate RHS term on each mesh
    mesh = obj.meshUnion(m);
    [ E, G ] = obj.matEvaluateFlux( mesh, fphys{m} );
    [ dFlux ] = obj.matEvaluateNumericalFlux( mesh, fphys{m}, obj.fext{m} );
    
    for i = 1:obj.Nvar
        [ obj.frhs{m}(:,:,i) ] = ...
            - mesh.rx.*( mesh.cell.Dr * E(:,:,i) ) ...
            - mesh.sx.*( mesh.cell.Ds * E(:,:,i) ) ...
            - mesh.ry.*( mesh.cell.Dr * G(:,:,i) ) ...
            - mesh.sy.*( mesh.cell.Ds * G(:,:,i) ) ...
            + ( mesh.cell.LIFT * ( mesh.Js .* dFlux(:,:,i) ))./ mesh.J;
    end
end

obj.matEvaluateSourceTerm( fphys );

end