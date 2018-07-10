%> @brief Evaluating the RHS term for the 2d problem
%> @details The function calculate the RHS term on each mesh
function matEvaluateRHS( obj, fphys )

obj.advectionSolver.evaluateAdvectionRHS( fphys );
obj.viscositySolver.matEvaluateRHS( fphys );
obj.matEvaluateSourceTerm( fphys );

end