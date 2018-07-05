classdef LSWEAbstract3d < handle
    
    properties ( Constant )
        %> num of 2d physical field
        Nfield2d = 5 % [ eta, HU, HV, H, Z ]
        %> num of 3d physical field
        Nfield3d = 7 % [ u, v, w, tau_x, tau_y, H, eta ]
        %> min water depth
        hmin = 1e-4
        %> gravity acceleration
        gra = 9.81
    end
    
    properties ( Abstract, Constant )
        startTime, finalTime
        %> output time interval
        outputTimeInterval
        %> case name
        casename
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
        %> viscosity
        miu
        %> linear slip parameter
        K
        %> output file object
        outputFile
    end
    
    properties ( SetAccess = public )
        %> time interval
        dt
    end
    
    methods ( Access = public )
        function obj = LSWEAbstract3d
        end% func
        
        function matSolve( obj )
            matEvaluateRK45( obj );
        end
        
        initPhysFromOptions( obj, mesh2d, mesh3d );
        drawVerticalSlice( obj, cellId_2d, nodeId_2d, field3d )
        AnimationSurfaceLevel( obj );
    end
    
    methods ( Access = protected )
        matUpdateExternalField( obj, time, fphys2d, fphys3d );
        matUpdateOutputResult( obj, time, fphys2d, fphys3d );
        [ fphys3d ] = matEvaluateVerticalVelocity( obj, mesh3d, fphys2d, fphys3d );
    end
    
    methods ( Abstract, Access = protected )
        %> evaluate boundary flux for 2d PCE
        frhs2d = EvaluatePCE2d_BoundaryKernel( obj, edge, fphys2d, fext2d )
        %> evaluate horizontal flux for 3d momentum equations
        frhs3d = Evaluate3d_HorizontalBoundaryKernel( obj, edge, fphys3d, fext3d )
    end
end