function matEvaluateSourceTerm( obj, fphys )

    matEvaluateSourceTerm@SWEAbstractCB2d( obj, fphys );
    
    % frhs = frhs + WindTerm    
    obj.windSolver.evaluateWindTermRHS(obj,fphys);
    
end
