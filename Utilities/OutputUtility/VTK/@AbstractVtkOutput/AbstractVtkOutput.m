classdef AbstractVtkOutput < AbstractOutputFile
    properties (SetAccess = protected)
        %> local node connection
        SEToV
        %> num of subcell
        Ncell
        %> num of points in each cell
        Np
        %> 
        CellVertList
        %> num of connection
        Ncon
        %> points coordinate
        Points
        %> num of points
        Npoint
        %> cell type
        ctype
    end

    methods ( Access = public )
        function obj = AbstractVtkOutput( casename, Nfield, dt )
            obj = obj@AbstractOutputFile( casename, Nfield, dt );
        end
        
        outputResult( obj, time, field );
        %> output result at specific step
        outputStepResult( obj, step, field );
    end
end