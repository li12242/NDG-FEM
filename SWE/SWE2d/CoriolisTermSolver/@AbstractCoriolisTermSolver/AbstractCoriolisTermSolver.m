classdef AbstractCoriolisTermSolver < handle
    
    properties(Constant)
    end
        
    methods(Abstract)
        frhs = evaluateCoriolisTermRHS( physClass, fphys )
    end
    
end

