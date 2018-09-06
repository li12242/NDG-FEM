classdef VtkOutput1d < AbstractVtkOutput
    %VTKOUTPUT1D Summary of this class goes here
    %   Detailed explanation goes here
    
    methods (Access = public)
        function obj = VtkOutput1d( casename, Nfield, dt )
            obj = obj@AbstractVtkOutput( casename, Nfield, dt );
        end

        function initFromMesh( obj, mesh )
            if ~isdir(obj.casename)
                mkdir(obj.casename);
            end
        end

        function closeOutputFile( obj )
        end

        % read output
        function readOutputResult( obj, timeStep )
        end
    end
end

