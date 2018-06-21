classdef VtkOutput2d < AbstractVtkOutput
    
methods (Access = public)
    function obj = VtkOutput2d( casename, Nfield, dt )
        obj = obj@AbstractVtkOutput( casename, Nfield, dt );
    end

    initFromMesh( obj, mesh );
    
    function closeOutputFile( obj )
    end
    
    % read output
    readOutputResult( obj, timeStep )
end
end