classdef WW1 < SDBAbstractTest

    properties( Constant )
        hmin = 1e-4
        gra = 9.81
    end
    
    properties( Constant )
        xlim = [-1, 1]
        ylim = [-1, 1]
    end
    
    properties( Constant )
        zm = 0
        zp = 1
        hm = 0.5
        hp = 0.5
    end
    
    methods
        function obj = WW1()
            obj = obj@SDBAbstractTest();
            mesh = obj.makeUniformMesh( obj.xlim, obj.ylim );
            obj.initPhysFromOptions( mesh );
        end        
    end
end
