classdef NdgNcOutputFileByStep < NdgNcOutputFile
    
    properties( SetAccess = protected )
        outputStepInterval
    end
    
    methods
        function obj = NdgNcOutputFileByStep( filename, ncdim, ncvar, outputStepInterval )
            obj = obj@NdgNcOutputFile( filename, ncdim, ncvar );
            obj.outputStepInterval = outputStepInterval;
        end
        
        function outputVar(obj, computeTime, varIndex, varargin)
            if ( obj.computeStep > (obj.outputStep * obj.outputStepInterval) )
                outputVar@NdgNcOutputFile(obj, varIndex, varargin{:});
            end
            obj.computeStep = obj.computeStep + 1;
        end
        
    end
    
end

