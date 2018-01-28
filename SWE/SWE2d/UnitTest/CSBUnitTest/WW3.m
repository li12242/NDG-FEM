classdef WW3 < CSBAbstractTest
    %WW3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties( Constant )
        h0 = 0.5;
        xlim = [ 0, 1];
        ylim = [-1, 1];
        hmin = 1e-4;
        gra = 9.81;
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
        
        function testPostFunc( obj, ntimes )
            for n = 1:ntimes
                %fphys = obj.matEvaluateLimiter( obj.fphys );
                [ fphys ] = obj.matEvaluatePostFunc( obj.fphys );
                obj.meshUnion.draw( fphys{1}(:,:,1) );
                drawnow;
            end
        end
        
        function testLimiter( obj, ntimes )
            for n = 1:ntimes
                fphys = obj.matEvaluateLimiter( obj.fphys );
                obj.meshUnion.draw( fphys{1}(:,:,1) );
                drawnow;
            end
        end
    end
    
    methods( Access = protected )
        function fphys = setInitialField( obj )
            fphys = cell(1, 1);
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                fphys{1} = zeros( mesh.cell.Np, mesh.K, obj.Nfield );
                fphys{1}(:,:,4) = obj.bx * mesh.x + obj.by * mesh.y;
                fphys{1}(:,:,1) = max( obj.h0 - fphys{1}(:,:,4), 0 );
            end
        end
    end
end

function [ mesh ] = makeUniformMesh( test )
bctype = [...
    NdgEdgeType.ZeroGrad, ...
    NdgEdgeType.ZeroGrad, ...
    NdgEdgeType.ZeroGrad, ...
    NdgEdgeType.ZeroGrad];

mesh = makeUniformTriMesh(1, test.xlim, test.ylim, 1, 10, bctype);
end% func


