classdef NdgSolver < handle
    
    properties( SetAccess = protected )
        %> enumeration for determining slope limiter (none or other types)
        limiterType
    end
    
    properties( SetAccess = protected )
        %> temporal discrete type (Euler, RK45 et al.)
        temporalDiscreteType
        %> enumeration for determining time interval-dt
        timeIntervalType
        %> time interval
        timeInterval
    end
    
    properties(SetAccess = protected)
        %> enumeration for open boundary types (function or NetCDF files)
        obcType
        %> enumeration types for updating boundary conditions
        obcIntervalType
        %> time interval for updating boundary conditions
        obcTimeInterval
        %> time step for updating boundary conditions
        obcStepInterval
        %> open boundary files
        obcFileName
    end
    
    properties(SetAccess = protected)
        %> name of output NetCDF file
        outputNetcdfCaseName
        %> enumeration types for output results
        outputIntervalType
        %> output step interval
        outputInterval
    end
    
    methods
        
        function obj = NdgSolver()
        end% func
        
        %> check whether all the properties in the class are assigned.
        checkPhysProperties( obj )
    end
    
end

