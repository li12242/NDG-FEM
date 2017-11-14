function [ outputfile ] = getOutputFile( FileType, outputIntervalType, varargin )

switch FileType
    case NdgIOFileType.None
        
    case NdgIOFileType.NetCDF
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
    case NdgIOIntervalType.DeltaTime
        outputfile = NdgNcOutputFileByTime( varargin{:} );
    case NdgIOIntervalType.DeltaStep
        outputfile = NdgNcOutputFileByStep( varargin{:} );
    otherwise
        msgID = [mfilename, ':getOutputIntervalTypeInvalid'];
        msgtext = 'The output interval type is unknown.';
        ME = MException(msgID, msgtext);
        throw(ME);
end

end
