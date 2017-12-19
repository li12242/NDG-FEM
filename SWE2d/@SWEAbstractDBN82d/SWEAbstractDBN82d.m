classdef SWEAbstractDBN82d < SWEAbstractDB2d

    
    properties(Constant)
        %> Physical field - {h, hu, hv, b, bx, by, Uwind, Vwind}
        Nfield = 8
    end
    
    properties
        windSolver
    end
    
    methods
        function obj = SWEAbstractDBN82d()
            obj = obj@SWEAbstractDB2d();
        end
        
        initPhysFromOptions( obj, mesh );
    end
    
    methods( Access = protected )
       matEvaluateSourceTerm( obj, mesh, fphys ) 
    end
    
end


