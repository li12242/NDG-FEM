function matEvaluateSourceTerm( obj, fphys )
% frhs = frhs + BottomTerm
obj.matEvaluateTopographySourceTerm( fphys );
% frhs = frhs + CoriolisTerm
obj.coriolisSolver.evaluateCoriolisTermRHS(obj, fphys);
% frhs = frhs + FrictionTerm
obj.frictionSolver.evaluateFrictionTermRHS(obj, fphys);
% frhs = frhs + WindTerm
obj.windSolver.evaluateWindTermRHS(obj, fphys);
end