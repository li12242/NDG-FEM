classdef NcOutput < AbstractOutputFile

    properties ( SetAccess = protected )
        %> output NetCDF file
        ncfile
        timeVarableId
        fieldVarableId
        filename
        vtkOutput
    end
    
    methods
        function obj = NcOutput( casename, Nfield, dt ) 
            obj = obj@AbstractOutputFile( casename, Nfield, dt );
        end

        %> create NetCDF output file
        initFromMesh( obj, mesh );
        %> output result
        outputResult( obj, time, field );

        [ field ] = readOutputResult( obj, step );

        writeResultToVtk( obj, step, field );
        writeOutputResultToVtk( obj, step );

        function closeOutputFile( obj )
            obj.ncfile.delete();
            obj.outputTime = ncread( obj.filename, 'time' );
        end
    end
    
end

