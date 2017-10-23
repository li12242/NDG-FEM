classdef NdgPhys < handle
    
    properties(SetAccess = protected)
        %> number of physical field
        Nfield
        %> mesh objects
        mesh
    end
    
    properties
        %> field values, dimensions - [Np, K, Nfield]
        fval
        %> final time
        finalTime
    end
    
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
        outputNetcdfFileName
        %> enumeration types for output results
        outputIntervalType
        %> output step interval
        outputInterval
    end
    
    methods
        function obj = NdgPhys(Nfield, mesh)
            obj.Nfield = Nfield;
            obj.mesh = mesh;
        end% func
        
        %> solve the problem with given condition
        function mxSolve( obj, mexSolveFunc )
            checkPhysProperties( obj );
            mexSolveFunc( obj );
        end
        
        %> check whether all the properties in the class are assigned.
        checkPhysProperties( obj )
    end
    
    methods(Access = protected)
               
    end
    
end

