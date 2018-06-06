classdef TidalRiverChannel2d < SWEConventional2d
    
    properties(Constant)
        %> water depth threshold
        hmin = 1e-4; 
        %> gravity acceleration
        gra = 9.8; 
        
        %> channel depth
        H = 320;
        %> wave period
        T = 500;
        %> delta = amplitude / H
        delta = 0.01
        
        %> domain width
        ChWidth = 500;
        %> river discharge rate
        dischargeRate = - 320 * sqrt( 9.8/320 );
    end
    
    properties( Abstract, Constant )
        %> left position of channel
        ChLeft
        %> right position of channel
        ChRight
    end
    
    methods( Access = public )
        function obj = TidalRiverChannel2d( N, M )
            obj = obj@SWEConventional2d();
            mesh = obj.makeUniformMesh( N, M );
            obj.initPhysFromOptions( mesh );
        end
    end
    
    methods( Access = protected, Static  )
        %> set open boundary condition
        obtype = setOpenBoundaryCondition( )
    end
    
    methods( Access = private )
        function mesh = makeUniformMesh( obj, N, M )
            obtype = obj.setOpenBoundaryCondition();
            bctype = [ NdgEdgeType.SlipWall, NdgEdgeType.SlipWall, ...
                obtype(1), obtype(2) ];
            mesh = makeUniformQuadMesh(N, ...
                [ obj.ChLeft, obj.ChRight ], [0, obj.ChWidth], M, 1, bctype);
        end
    end
    
    methods( Access = protected )
        function fphys = setInitialField( obj )
            fphys = cell( obj.Nmesh, 1 );
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                fphys{m} = zeros( mesh.cell.Np, mesh.K, obj.Nfield );
                fphys{m}(:, :, 1) = + obj.H;
                fphys{m}(:, :, 4) = - obj.H;
            end
        end
        
        function matUpdateExternalField( obj, time, fphys )
            for m = 1:obj.Nmesh
                a = 1 - 1./exp( time/obj.T/2 );
                mesh = obj.meshUnion(m);
                [ obj.fext{m}(:,:,1), obj.fext{m}(:,:,2) ] ...
                    = obj.setOBC( mesh.x, time );
                
                if time < obj.T/2
                    obj.fext{m}(:, :, 1) = a.* ( obj.fext{m}(:, :, 1) - obj.H ) + obj.H;
                    obj.fext{m}(:, :, 2) = a.* obj.fext{m}(:, :, 2);
                end
            end
        end
        
        %> Set open boundary condition
        function [h, hu] = setOBC( obj, x, time )
            h = zeros( size(x) );
            hu = zeros( size(x) );
            
            w = 2 * pi / obj.T;
            c = sqrt( obj.gra * obj.H );
            k = w/c;
            
            temp = cos( k .* x(:, 1) - w * time );
            eta = obj.delta * obj.H;
            h(:, 1) = eta * temp + obj.H;
            u = eta * sqrt( obj.gra/obj.H ) * temp;
            hu(:, 1) = h(:, 1) .* u;
            
            hu(:, end) = obj.dischargeRate;
        end
    end
end

