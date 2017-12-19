classdef AbstractFrictionTermSolver < handle
    
    properties(Constant)

    end
        
    methods(Abstract)
        frhs = evaluateFrictionTermRHS( physClass, fphys )
    end
    
end

