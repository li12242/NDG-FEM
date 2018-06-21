classdef VtkOutput3d < VtkOutput2d
    methods (Access = public)
        function obj = VtkOutput3d( casename, Nfield, dt )
            obj = obj@VtkOutput2d( casename, Nfield, dt );
        end

        initFromMesh( obj, mesh );
        % drawResult( obj )
        
        readOutputResult( obj, timeStep )
    end
end