classdef WDTest1 < WDTestAbstract
    %WDTEST1 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties( Constant )
        hmin = 1e-4
        gra = 9.81
    end
    
    methods
        function obj = WDTest1(N)
            obj = obj@WDTestAbstract();
            [ mesh ] = makeUniformMesh( N );
            obj.initPhysFromOptions( mesh );
        end
    end
    
    methods( Access = protected )
        function fphys = setInitialField( obj )
            fphys = cell( obj.Nmesh, 1 );
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                fphys{m} = zeros(mesh.cell.Np, mesh.K, obj.Nfield);
                fphys{m}(:, 1, 1) = 1;
                fphys{m}(:, 2, 1) = 1 - mesh.x(:, 2);
            end
        end% func
    end
    
end

function [ mesh ] = makeUniformMesh( N )
xlim = [-1, 1];
bcType = [NdgEdgeType.ZeroGrad, NdgEdgeType.ZeroGrad];
[ mesh ] = makeUniformMesh1d( N, xlim, 2, bcType );
end