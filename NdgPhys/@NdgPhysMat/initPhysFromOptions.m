%> @brief Initialize the soler with the options, fphysical field and solver methods.
function initPhysFromOptions( obj, mesh )
    % call the superclass methods
    initPhysFromOptions@NdgPhys( obj, mesh );
    
    % set the physical field for the NdgPhysMat solver
    for m = 1:obj.Nmesh
        Np = obj.meshUnion(m).cell.Np;
        K = obj.meshUnion(m).K;
        obj.frhs{m} = zeros( Np, K, obj.Nvar );
        if ~isempty( mesh(m).BoundaryEdge )
            Nfp = mesh(m).BoundaryEdge.Nfp;
            Ne = mesh(m).BoundaryEdge.Ne;
            obj.fext{m} = zeros( Nfp, Ne, obj.Nfield );
        end
    end

    % Setup the output NetCDF file object
    obj.outputFile = initOutput( obj, mesh );
    
    % initilize the solver methods
    [ obj.advectionSolver, obj.viscositySolver ] = initSolver( obj );

    %> set the final time
    obj.ftime = obj.getOption('finalTime');
    if obj.option.isKey('CFL')
        obj.cfl = obj.getOption('CFL');
    elseif obj.option.isKey('cfl')
        obj.cfl = obj.getOption('cfl');
    elseif obj.option.isKey('Cfl')
        obj.cfl = obj.getOption('Cfl');
    else
        obj.cfl = 1;
    end
    %> choose the limiter
    [ obj.limiter ] = initSlopeLimiter( obj );
end% func