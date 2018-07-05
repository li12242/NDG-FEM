function initPhysFromOptions( obj, mesh2d, mesh3d )

    % set mesh object
    obj.mesh2d = mesh2d;
    obj.mesh3d = mesh3d;
    obj.Nmesh = numel( mesh3d );

    % set the physical field for the NdgPhysMat solver
    for m = 1:obj.Nmesh
        Np = obj.mesh3d(m).cell.Np;
        K = obj.mesh3d(m).K;
        obj.frhs3d{m} = zeros( Np, K, 2 );
        if ~isempty( mesh3d(m).BoundaryEdge )
            Nfp = mesh3d(m).BoundaryEdge.Nfp;
            Ne = mesh3d(m).BoundaryEdge.Ne;
            obj.fext3d{m} = zeros( Nfp, Ne, obj.Nfield3d );
        end
        
        Np = obj.mesh2d(m).cell.Np;
        K = obj.mesh2d(m).K;
        obj.frhs2d{m} = zeros( Np, K, 1 );
        if ~isempty( mesh2d(m).BoundaryEdge )
            Nfp = mesh2d(m).BoundaryEdge.Nfp;
            Ne = mesh2d(m).BoundaryEdge.Ne;
            obj.fext2d{m} = zeros( Nfp, Ne, obj.Nfield2d );
        end
    end

    % Setup the output NetCDF file object
    initOutput( obj, mesh2d, mesh3d );
end