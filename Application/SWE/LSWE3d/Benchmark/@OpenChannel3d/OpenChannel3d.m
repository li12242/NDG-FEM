classdef OpenChannel3d < LSWEAbstract3d
    %OPENCHANNEL3D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Constant)
        ChLength = 28e3 * 0.8;
        ChWidth = 500;
        
        %> channel depth
        H = 320;
        startTime = 0;
        finalTime = 5000;
        dt = 1;
        eta = 0.5;
        
        %> viscosity
        miu = 0.0015;
        K = 0.01;
    end
    
    methods
        function obj = OpenChannel3d( N, Nz, M, Mz )
            [ obj.mesh2d, obj.mesh3d ] = makeChannelMesh( obj, N, Nz, M, Mz );
            obj.initPhysFromOptions( obj.mesh2d, obj.mesh3d );
        end
    end
    
    methods ( Access = protected )
        function [fphys2d, fphys3d] = setInitialField( obj )
            fphys2d = cell( obj.Nmesh, 1 );
            fphys3d = cell( obj.Nmesh, 1 );
            for m = 1 : obj.Nmesh
                mesh2d = obj.mesh2d(m);
                mesh3d = obj.mesh3d(m);
                fphys2d{m} = zeros( mesh2d.cell.Np, mesh2d.K, obj.Nfield2d );
                fphys3d{m} = zeros( mesh3d.cell.Np, mesh3d.K, obj.Nfield3d );
                
                % surface elevation
                % fphys2d{m}(:, :, 1) = obj.eta * ( mesh2d.x ./ obj.ChLength + 0.5);
                fphys2d{m}(:, :, 1) = obj.eta * sin( ( mesh2d.x ./ obj.ChLength + 0.5) * pi );
                % bottom elevation
                fphys2d{m}(:, :, 5) = - obj.H;
                % water depth
                fphys2d{m}(:, :, 4) = fphys2d{m}(:, :, 1) - fphys2d{m}(:, :, 5);
                
                % water depth
                fphys3d{m}(:, :, 5) = mesh3d.Extend2dField( fphys2d{m}(:, :, 4) );
                % bottom elevation
                fphys3d{m}(:, :, 6) = mesh3d.Extend2dField( fphys2d{m}(:, :, 1) );
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
mesh3d.BottomEdge = NdgBottomEdge3d( mesh3d, 1 );
mesh3d.BoundaryEdge = NdgHaloEdge3d( mesh3d, 1 );

end
