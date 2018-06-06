%> @brief Initialize the soler with the options, fphysical field and solver methods.
function initPhysFromOptions( obj, mesh )
    % call the superclass methods
    initPhysFromOptions@NdgPhys( obj, mesh );
    
    % set the physical field for the NdgPhysMat solver
    for m = 1:obj.Nmesh
        Np = obj.meshUnion(m).cell.Np;
        K = obj.meshUnion(m).K;
        obj.frhs{m} = zeros( Np, K, obj.Nvar );
%         Nfp = mesh(m).BoundaryEdge.Nfp;
%         Ne = mesh(m).BoundaryEdge.Ne;
%         obj.fext{m} = zeros( Nfp, Ne, obj.Nfield );
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

% function [ ncfile ] = makeBasicOutputNetcdfFile( phys, ind, mesh )
%     dimTime = NdgNcDim('Nt', 0);
%     dimK = NdgNcDim('K', mesh.K);
%     dimNp = NdgNcDim('Np', mesh.cell.Np);
%     dimNfield = NdgNcDim('Nvar', phys.Nvar);
% 
%     varTime = NdgNcVar('time', dimTime, enumNcData.NC_DOUBLE );
%     varField = NdgNcVar('fphys', [dimNp, dimK, dimNfield, dimTime], enumNcData.NC_DOUBLE);
% 
%     filename = [ phys.getOption('outputNetcdfCaseName') , '.', num2str(ind),...
%         '-', num2str( phys.Nmesh ), '.nc' ];
% 
%     intervalType = phys.getOption('outputIntervalType');
%     switch intervalType
%         case enumOutputInterval.DeltaTime
%             interval = phys.getOption('outputTimeInterval');
%         case enumOutputInterval.DeltaStep
%             interval = phys.getOption('outputStepInterval');
%         otherwise
%             msgID = [mfilename, ':outputIntervalTypeInvalid'];
%             msgtext = 'The output type is unknow.';
%             throw( MException(msgID, msgtext) );
%     end
%     ncfile = getOutputFile( enumOutputFile.NetCDF, intervalType, filename, ...
%         [dimTime, dimK, dimNp, dimNfield], [varTime, varField], interval);
% end% func