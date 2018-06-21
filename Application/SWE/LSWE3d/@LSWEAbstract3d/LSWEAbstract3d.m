classdef LSWEAbstract3d < handle
    
    properties ( Constant )
        %> num of 2d physical field
        Nfield2d = 5 % [ eta, HU, HV, H, Z ]
        %> num of 3d physical field
        Nfield3d = 6 % [ u, v, tau_x, tau_y, H, eta ]
        %> min water depth
        hmin = 1e-4
        %> gravity acceleration
        gra = 9.81
        
        
    end
    
    properties ( Abstract, Constant )
        startTime, finalTime
        %> time interval
        dt
        %> viscosity
        miu;
        %> slip friction parameter
        K
    end
    
    properties ( SetAccess = protected )
        %> num of mesh
        Nmesh
        %> physical field
        fphys2d, frhs2d, fext2d
        %> physical field
        fphys3d, frhs3d, fext3d
        %> horizontal mesh
        mesh2d
        %> vertical extended mesh
        mesh3d
    end

    properties ( SetAccess = protected )
        %> output file object
        outputFile
    end

    methods ( Access = public )
        function obj = LSWEAbstract3d
        end% func
        
        function matSolve( obj )
            matEvaluateRK45( obj );
        end
        
        initPhysFromOptions( obj, mesh2d, mesh3d );
    end
    
    methods ( Access = protected )
        matUpdateOutputResult( obj, time, fphys )
    end
    
    methods ( Abstract, Access = protected )
        % set initilize physical field
        fphys = setInitialField( obj )
    end
end