function [ outputfile ] = getOutputFile( FileType, outputIntervalType, varargin )

switch FileType
    case enumOutputFile.None
        
    case enumOutputFile.NetCDF
        outputfile = getOutputNcFile( outputIntervalType, varargin{:} );
        
        
    otherwise
        msgID = [mfilename, ':getOutputFileTypeInvalid'];
        msgtext = 'The output file type is unknown.';
        ME = MException(msgID, msgtext);
        throw(ME);
end

end

function outputfile = getOutputNcFile( outputIntervalType, varargin )

switch outputIntervalType
    case enumOutputInterval.DeltaTime
        outputfile = NdgNcOutputFileByTime( varargin{:} );
    case enumOutputInterval.DeltaStep
        outputfile = NdgNcOutputFileByStep( varargin{:} );
    otherwise
        msgID = [mfilename, ':getOutputIntervalTypeInvalid'];
        msgtext = 'The output interval type is unknown.';
        ME = MException(msgID, msgtext);
        throw(ME);
end

end
