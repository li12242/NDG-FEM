classdef NdgAbstractOutputFile < handle
    
    properties
    end
    
    methods
        function obj = NdgAbstractOutputFile()
        end
    end
    
    methods( Abstract )
        %> output the field into files
        outputVar( obj, computeTime, varIndex, varargin)
        %> output the final result
        outputFinalVar( obj, varIndex, varargin )
        %> initialize all the properties of the output file
        initOutputFile( obj )
    end
    
end

