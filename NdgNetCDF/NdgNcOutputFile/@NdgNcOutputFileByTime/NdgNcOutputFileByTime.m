classdef NdgNcOutputFileByTime < NdgNcOutputFile
    
    properties( SetAccess = protected )
        outputTimeInterval
    end
    
    methods
        function obj = NdgNcOutputFileByTime( filename, ncdim, ncvar, outputDt )
            obj = obj@NdgNcOutputFile( filename, ncdim, ncvar );
            obj.outputTimeInterval = outputDt;
        end
        
        function outputVar(obj, computeTime, varIndex, varargin)
            if ( computeTime > (obj.outputStep * obj.outputTimeInterval) )
                outputVar@NdgNcOutputFile(obj, varIndex, varargin{:});
            end
            obj.computeStep = obj.computeStep + 1;
        end
        
    end
    
end
