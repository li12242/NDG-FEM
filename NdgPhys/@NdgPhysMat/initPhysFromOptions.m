%> @brief Initialize the soler with the options, fphysical field and solver methods.
function initPhysFromOptions( obj, mesh )
    % call the superclass methods
    initPhysFromOptions@NdgPhys( obj, mesh );
    % set the physical field for the NdgPhysMat solver
    for m = 1:obj.Nmesh
        Np = obj.meshUnion(m).cell.Np;
        K = obj.meshUnion(m).K;
        obj.frhs{m} = zeros( Np, K, obj.Nvar );
        obj.fext{m} = zeros( Np, K, obj.Nfield );
    end

    % Setup the output NetCDF file object
    obj.outputNcFile = [];
    for n = 1:obj.Nmesh
        [ ncfile ] = makeBasicOutputNetcdfFile( obj, n, obj.meshUnion(n) );
        obj.outputNcFile = [obj.outputNcFile, ncfile];
    end
    % initilize the solver methods
    [ obj.advectionSolver ] = setAdvectionSolver( obj );
    [ obj.viscositySolver ] = setViscositySolver( obj );

    %> set the final time
    obj.ftime = obj.getOption('finalTime');
    if obj.option.isKey('CFL')
        obj.cfl = obj.getOption('CFL');
    elseif obj.option.isKey('cfl')
        obj.cfl = obj.getOption('cfl');
    else
        obj.cfl = 1;
    end
    %> choose the limiter
    [ obj.limiter ] = setSlopeLimiter( obj );
end% func

%> Set the slope limiter object from option
function [ limiter ] = setSlopeLimiter( obj )
    if obj.option.isKey('limiterType')
        switch obj.getOption('limiterType')
            case NdgLimiterType.None
                limiter = NdgNonLimiter( obj.meshUnion );
            case NdgLimiterType.Vert
                limiterCreator = choseFuncByMeshDim( obj, ...
                    @NdgVertLimiter1d, ...
                    @NdgVertLimiter2d );
                limiter = limiterCreator( obj.meshUnion );
            case NdgLimiterType.TVB
                limiterCreator = choseFuncByMeshDim( obj, ...
                    @NdgTVB1d, ...
                    @NdgTVB2d );
                M = obj.getOption('limiterParameter');
                limiter = limiterCreator( obj.meshUnion, M );
            case NdgLimiterType.BJ
                limiterCreator = choseFuncByMeshDim( obj, ...
                    @NdgBJ1d, ...
                    @NdgBJ2d );
                limiter = limiterCreator( obj.meshUnion );
        end% switch
    else
        limiter = NdgNonLimiter( obj.meshUnion );
    end
end% func

function [ viscositySolver ] = setViscositySolver( obj )
    if obj.option.isKey('viscosityType')
        switch obj.getOption('viscosityType')
            case NdgViscosityType.None
                viscositySolver = NdgNonVisSolver( obj );
            case NdgViscosityType.LDG
        end% switch
    else
        viscositySolver = NdgNonVisSolver( obj );
    end
end% func

function [ advectionSolver ] = setAdvectionSolver( obj )
    integralType = obj.getOption('integralType');
    equationType = obj.getOption('equationType');
    
    switch integralType
        case NdgDiscreteIntegralType.QuadratureFree
            
            switch equationType
                case NdgDiscreteEquationType.Strong
                    advectionCreator = choseFuncByMeshDim( obj, ...
                        @NdgQuadFreeStrongFormAdvSolver1d, ...
                        @NdgQuadFreeStrongFormAdvSolver2d);
                case NdgDiscreteEquationType.Weak
                    advectionCreator = choseFuncByMeshDim( obj, ...
                        @NdgQuadFreeWeakFormAdvSolver1d, ...
                        @NdgQuadFreeWeakFormAdvSolver2d);
            end
            
        case NdgDiscreteIntegralType.GaussQuadrature
            
            switch equationType
                case NdgDiscreteEquationType.Strong
                    advectionCreator = choseFuncByMeshDim( obj, ...
                        @NdgGaussQuadStrongFormAdvSolver1d, ...
                        @NdgGaussQuadStrongFormAdvSolver2d);
                case NdgDiscreteEquationType.Weak
                    advectionCreator = choseFuncByMeshDim( obj, ...
                        @NdgGaussQuadWeakFormAdvSolver1d, ...
                        @NdgGaussQuadWeakFormAdvSolver2d);
            end
            
        otherwise
            msgID = [mfilename, ':inputIntegralTypeInvalid'];
            msgtext = ['please set the "integralType" in the option ', ...
                'from the existing NdgDiscreteIntegralType class names.'];
            throw( MException(msgID, msgtext) );
    end% switch
    
    advectionSolver = advectionCreator( obj );
end% func

function [ func ] = choseFuncByMeshDim( phys, func1d, func2d )
    switch phys.meshUnion(1).type
        case NdgMeshType.OneDim
            func = func1d;
        case NdgMeshType.TwoDim
            func = func2d;
    end
end% func

function [ ncfile ] = makeBasicOutputNetcdfFile( phys, ind, mesh )
    dimTime = NdgNcDim('Nt', 0);
    dimK = NdgNcDim('K', mesh.K);
    dimNp = NdgNcDim('Np', mesh.cell.Np);
    dimNfield = NdgNcDim('Nvar', phys.Nvar);

    varTime = NdgNcVar('time', dimTime, NdgNcDataType.NC_DOUBLE );
    varField = NdgNcVar('fphys', [dimNp, dimK, dimNfield, dimTime], NdgNcDataType.NC_DOUBLE);

    filename = [ phys.getOption('outputNetcdfCaseName') , '.', num2str(ind),...
        '-', num2str( phys.Nmesh ), '.nc' ];

    intervalType = phys.getOption('outputIntervalType');
    switch intervalType
        case NdgIOIntervalType.DeltaTime
            interval = phys.getOption('outputTimeInterval');
        case NdgIOIntervalType.DeltaStep
            interval = phys.getOption('outputStepInterval');
        otherwise
            msgID = [mfilename, ':outputIntervalTypeInvalid'];
            msgtext = 'The output type is unknow.';
            throw( MException(msgID, msgtext) );
    end
    ncfile = getOutputFile( NdgIOFileType.NetCDF, intervalType, filename, ...
        [dimTime, dimK, dimNp, dimNfield], [varTime, varField], interval);
end% func