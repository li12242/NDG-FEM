classdef SWEAbstractDBN62d < SWEAbstractDB2d

    
    properties(Constant)
        %> Physical field - {h, hu, hv, b, bx, by}
        Nfield = 6
    end
    
    methods
        function obj = SWEAbstractDBN62d()
            obj = obj@SWEAbstractDB2d();
        end        
    end
    
end

