classdef NdgNcOutputFile < NdgNcFile
    %NDGOUTPUTFILE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        
        function obj = NdgNcOutputFile( varargin )
            obj = obj@NdgNcFile( varargin{:} );
        end
        
        function outputTotalVar(obj, varIndex, value)
            netcdf.putVar(obj.ncid, obj.ncVar(varIndex).id, value);
        end
        
        function outputVar(obj, varIndex, value, outputStep)
            Ndim = numel(obj.ncVar(varIndex).dims);
            startInd = [zeros(1, Ndim-1), outputStep];
            countInd = [obj.ncVar(varIndex).dims(1:end-1).length, 1];
            netcdf.putVar(obj.ncid, obj.ncVar(varIndex).id, ...
                startInd, countInd, value);
        end
    end
    
end