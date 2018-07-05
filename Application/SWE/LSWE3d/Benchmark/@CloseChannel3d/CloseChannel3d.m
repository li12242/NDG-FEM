classdef CloseChannel3d < LSWEAbstract3d
    %OPENCHANNEL3D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Constant)
        %> channel length
        ChLength = 40e3;
        %> channel width
        ChWidth = 8e3;
        %> channel depth
        H = 12;
        %> start time
        startTime = 0;
        %> final time
        finalTime = 24 * 3600;
        %> output time interval
        outputTimeInterval = 60;
        %> casename
        casename = 'FreeSurfaceSeichingWithHorizontalBottomClosedRectangularBasin';
        %> maximum amplitude
        eta = 0.25;
    end
    
    methods
        function obj = CloseChannel3d( N, Nz, M, Mz )
            % setup mesh domain
            [ obj.mesh2d, obj.mesh3d ] = makeChannelMesh( obj, N, Nz, M, Mz );
            % allocate boundary field with mesh obj 
            obj.initPhysFromOptions( obj.mesh2d, obj.mesh3d );
            % set initilize physical field
            [ obj.fphys2d, obj.fphys3d ] = obj.setInitialField( );
            % set vertical viscosity
            obj.miu = 0;
            %> linear slip parameter
            obj.K = 0;
            %> time interval
            obj.dt = 80;
        end 
        
        AnalysisResult2d( obj )
        AnalysisResult3d( obj )
    end
    
    methods ( Access = protected )
        %> evaluate boundary flux for 2d PCE
        frhs2d = EvaluatePCE2d_BoundaryKernel( obj, edge, fphys2d, fext2d )
        %> evaluate horizontal flux for 3d momentum equations
        frhs3d = Evaluate3d_HorizontalBoundaryKernel( obj, edge, fphys3d, fext3d )
        
        function matUpdateExternalField( obj, time, fphys2d, fphys3d )
        end
        
        %> set initial function
        function [fphys2d, fphys3d] = setInitialField( obj )
            fphys2d = cell( obj.Nmesh, 1 );
            fphys3d = cell( obj.Nmesh, 1 );
            for m = 1 : obj.Nmesh
                mesh2d = obj.mesh2d(m);
                mesh3d = obj.mesh3d(m);
                fphys2d{m} = zeros( mesh2d.cell.Np, mesh2d.K, obj.Nfield2d );
                fphys3d{m} = zeros( mesh3d.cell.Np, mesh3d.K, obj.Nfield3d );
                
                % surface elevation
                fphys2d{m}(:, :, 1) = obj.eta * cos( ( mesh2d.x ./ obj.ChLength ) * pi );
                % bottom elevation
                fphys2d{m}(:, :, 5) = - obj.H;
                % water depth
                fphys2d{m}(:, :, 4) = obj.H;
                
                % water depth
                fphys3d{m}(:, :, 6) = mesh3d.Extend2dField( fphys2d{m}(:, :, 4) );
                % bottom elevation
                fphys3d{m}(:, :, 7) = mesh3d.Extend2dField( fphys2d{m}(:, :, 1) );
            end
        end
    end
end

function [mesh2d, mesh3d] = makeChannelMesh( obj, N, Nz, M, Mz )

bctype = [ ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall ];

mesh2d = makeUniformTriMesh( N, ...
    [ -obj.ChLength, 0 ], [0, obj.ChWidth], M, 2, bctype);

cell = StdPrismTri( N, Nz );
zs = zeros(mesh2d.Nv, 1); zb = zs - 1;
mesh3d = NdgExtendMesh3d( cell, mesh2d, zs, zb, Mz );
mesh3d.InnerEdge = NdgSideEdge3d( mesh3d, 1 );
mesh3d.BottomEdge = NdgBottomInnerEdge3d( mesh3d, 1 );
mesh3d.BoundaryEdge = NdgHaloEdge3d( mesh3d, 1 );
mesh3d.BottomBoundaryEdge = NdgBottomHaloEdge3d( mesh3d, 1 );
mesh3d.SurfaceBoundaryEdge = NdgSurfaceHaloEdge3d( mesh3d, 1 );
end
