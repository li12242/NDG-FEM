classdef VtkOutput2d < AbstractVtkOutput
    
methods (Access = public)
    function obj = VtkOutput2d( casename, Nfield, dt )
        obj = obj@AbstractVtkOutput( casename, Nfield, dt );
    end

    initFromMesh( obj, mesh );
    %outputResult( obj, time, field );
    function closeOutputFile( obj )
    end
    % drawResult( obj )
    readOutputResult( obj, timeStep )
end
end