classdef AbstractWindTermSolver < handle
    
    properties(Constant)

    end
    
    methods(Abstract)
        frhs = evaluateWindTermRHS( physClass, fphys )
    end
    
end

