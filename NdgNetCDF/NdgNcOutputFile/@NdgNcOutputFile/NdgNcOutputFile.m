classdef NdgNcOutputFile < NdgNcFile & NdgAbstractOutputFile
    
    properties( SetAccess = protected )
        computeStep = 0;
        outputStep = 0;
        outputStepInterval
    end
    
    methods
        
        function obj = NdgNcOutputFile( filename, ncdim, ncvar )
            obj = obj@NdgNcFile( filename, ncdim, ncvar );
            % create the NetCDF file
            obj.defineIntoNetcdfFile(); 
        end
        
        function initOutputFile( obj )
            obj.computeStep = 0;
            obj.outputStep = 0;
        end
        
        function outputTotalVar(obj, varIndex, value)
            netcdf.putVar(obj.ncid, obj.ncVar(varIndex).id, value);
        end
        
        function outputVar(obj, varIndex, varargin)
            
            for n = 1:numel(varIndex)
                ind = varIndex(n);
                value = varargin{n};
                Ndim = numel(obj.ncVar(ind).dims);
                startInd = [zeros(1, Ndim-1), obj.outputStep];
                countInd = [obj.ncVar(ind).dims(1:end-1).length, 1];
                netcdf.putVar(obj.ncid, obj.ncVar(ind).id, startInd, countInd, value);
            end
            obj.outputStep = obj.outputStep + 1;
        end
        
        function outputFinalVar( obj, computeTime, varIndex, varargin )
            for n = 1:numel(varIndex)
                ind = varIndex(n);
                value = varargin{n};
                Ndim = numel(obj.ncVar(ind).dims);
                startInd = [zeros(1, Ndim-1), obj.outputStep];
                countInd = [obj.ncVar(ind).dims(1:end-1).length, 1];
                netcdf.putVar(obj.ncid, obj.ncVar(ind).id, startInd, countInd, value);
            end
            obj.outputStep = obj.outputStep + 1;
        end
        
    end
    
end