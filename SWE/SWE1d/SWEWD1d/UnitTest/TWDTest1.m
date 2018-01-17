classdef TWDTest1 < TWDTestAbstract
    
    properties( Constant )
        hmin = 1e-4;
        gra = 9.81;
    end
    
    properties
        h0 = 10;
    end
    
    methods
        function obj = TWDTest1( N )
            obj = obj@TWDTestAbstract();
            [ mesh ] = makeUniformMesh( N, 2 );
            obj.initPhysFromOptions( mesh );
        end
    end
    
    methods( Access = protected )
        function fphys = setInitialField( obj )
            fphys = cell( obj.Nmesh, 1 );
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                fphys{1} = zeros(mesh.cell.Np, mesh.K, obj.Nfield);
                fphys{1}(:,1,1) = obj.h0;
            end
        end
    end
    
end

function [ mesh ] = makeUniformMesh( N, M )
xlim = [-1, 1];
bcType = [NdgEdgeType.ZeroGrad, NdgEdgeType.ZeroGrad];
[ mesh ] = makeUniformMesh1d( N, xlim, M, bcType );
end
