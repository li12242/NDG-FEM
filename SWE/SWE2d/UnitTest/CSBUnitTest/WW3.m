classdef WW3 < CSBAbstractTest
    %WW3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties( Constant )
        bx = 1;
        by = 0;
    end
    
    methods
        function obj = WW3()
            obj = obj@CSBAbstractTest();
            mesh = makeUniformMesh( obj );
            obj.initPhysFromOptions( mesh );
            obj.matEvaluatePostFunc( obj.fphys );
        end
    end
    
end

function [ mesh ] = makeUniformMesh( test )
bctype = [...
    NdgEdgeType.ZeroGrad, ...
    NdgEdgeType.ZeroGrad, ...
    NdgEdgeType.ZeroGrad, ...
    NdgEdgeType.ZeroGrad];

mesh = makeUniformQuadMesh(1, test.xlim, test.ylim, 1, 1, bctype);
end% func


