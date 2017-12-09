function matEvaluateSourceTerm( obj, fphys )

    for m = 1:obj.Nmesh
        mesh = obj.meshUnion(m);
        obj.frhs{m} = obj.frhs{m} + mxEvaluateSourceTopography2d...
            ( obj.gra, mesh.EToR, fphys{m} );
    end

    % frhs = frhs + CoriolisTerm
    obj.coriolisSolver.evaluateCoriolisTermRHS(obj,fphys);
    
    % frhs = frhs + FrictionTerm
    obj.frictionSolver.evaluateFrictionTermRHS(obj,fphys);
    

end
