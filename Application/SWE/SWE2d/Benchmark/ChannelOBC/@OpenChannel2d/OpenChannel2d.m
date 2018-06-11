classdef OpenChannel2d < SWEConventional2d
    
    properties(Constant)
        %> water depth threshold
        hmin = 1e-4; 
        %> gravity acceleration
        gra = 9.8; 
        %> channel depth
        H = 320;
        %> delta = amplitude / depth
        delta = 0.04;
        %> wave period
        T = 500;
        
        %> domain length = wave length * lambda
        ChLength = 28e3 * 0.8;
        %> domain width
        ChWidth = 500;
        %> mesh rotation angle
        theta = 0;
        %> mesh center
        xc = 0;
        %> mesh center
        yc = 0;
    end
    
    methods
        function obj = OpenChannel2d( N, M )
            obj = obj@SWEConventional2d();
            mesh = obj.makeUniformMesh( N, M );
            obj.initPhysFromOptions( mesh );
        end
        
        %> draw the numerical results at gauge points
        drawGaugePoints( obj, xg );

        %> show the paramters for determing the test problem
        showParameter( obj )
        
        %> check mechanical energy, excess mass, energy flux
        checkResult( obj )
    end
    
    methods( Access = protected, Static  )
        %> set open boundary condition
        obtype = setOpenBoundaryCondition( )
    end
    
    methods( Abstract, Access = protected )
        [h, hu] = setExactSolution( obj, x, time );
    end
    
    methods( Access = private, Sealed )
        function mesh = makeUniformMesh( obj, N, M )
            obtype = obj.setOpenBoundaryCondition();
            bctype = [ ...
                enumBoundaryCondition.SlipWall, ...
                enumBoundaryCondition.SlipWall, ...
                obtype(1), obtype(2)];
            
            mesh = makeUniformQuadMesh(N, ...
                [-obj.ChLength, 0], [0, obj.ChWidth], M, 1, bctype);
        end
    end
    
    methods( Access = protected )
                
        function fphys = setInitialField( obj )
            fphys = cell( obj.Nmesh, 1 );
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                fphys{m} = zeros( mesh.cell.Np, mesh.K, obj.Nfield );
%                 [ fphys{m}(:,:,1), fphys{m}(:,:,2) ] ...
%                     = obj.setExactSolution( mesh.x, 0.0 );
                fphys{m}(:, :, 1) = + obj.H;
                fphys{m}(:, :, 4) = - obj.H;
            end
        end
        
        function matUpdateExternalField( obj, time, fphys )
            for m = 1:obj.Nmesh
                edge = obj.meshUnion(m).BoundaryEdge;
                a = 1 - 1./exp( time/obj.T/2 );
                [ obj.fext{m}(:,:,1), obj.fext{m}(:,:,2) ] ...
                    = obj.setOBC( edge.xb, time );

                if time < obj.T/2
                    obj.fext{m}(:, :, 1) = a.* ( obj.fext{m}(:, :, 1) - obj.H ) + obj.H;
                    obj.fext{m}(:, :, 2) = a.* obj.fext{m}(:, :, 2);
                end
            end
        end
        
        %> Set open boundary condition
        function [h, hu] = setOBC( obj, x, time )
            w = 2 * pi / obj.T;
            c = sqrt( obj.gra * obj.H );
            k = w/c;
            
            temp = cos( k .* x - w * time );
            eta = obj.delta * obj.H;
            h = eta * temp + obj.H;
            u = eta * sqrt( obj.gra/obj.H ) * temp;
            hu = h .* u;
        end
                
    end% methods
    
end

