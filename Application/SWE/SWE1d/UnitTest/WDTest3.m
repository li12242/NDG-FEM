classdef WDTest3 < WDTestAbstract
    %WDTEST3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties( Constant )
        hmin = 1e-4
        gra = 9.81
    end
    
    methods
        function obj = WDTest3(N)
            obj = obj@WDTestAbstract();
            [ mesh ] = makeUniformMesh( N );
            obj.initPhysFromOptions( mesh );
        end
    end
    
    methods( Access = protected )
        function fphys = setInitialField( obj )
            fphys = cell( obj.Nmesh, 1 );
            eta = 0.5;
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                fphys{1} = zeros(mesh.cell.Np, mesh.K, obj.Nfield);
                
                fphys{1}(:, :, 3) = -( mesh.x );
                fphys{1}(:, 1, 1) = eta;
                fphys{1}(:, 2, 1) = max( 0, eta - abs( mesh.x(:,2) ));
            end
        end% func
    end
end

function [ mesh ] = makeUniformMesh( N )
xlim = [-1, 1];
bcType = [NdgEdgeType.ZeroGrad, NdgEdgeType.ZeroGrad];
[ mesh ] = makeUniformMesh1d( N, xlim, 2, bcType );
end

