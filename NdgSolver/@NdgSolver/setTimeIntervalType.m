%======================================================================
%> @brief Brief description of the function
%>
%> More detailed description.
%>
%> @param arg1 First argument
%> @param arg2 Second argument
%>
%> @retval out1 return value for the first output variable
%> @retval out2 return value for the second output variable
%======================================================================
%> This function is part of the NDGOM software. 
%> @author li12242, Tianjin University, li12242@tju.edu.cn
%======================================================================
function obj = setTimeIntervalType(obj, type, varargin)

try
    obj.timeIntervalType = NdgIntervalType( type );
catch 
    msgID = 'Phys:setTimeIntervalType';
    msgtext = ['The input time discrete type should be an ',...
        'enumeration from NdgIntervalType.';];
    ME = MException(msgID, msgtext);
    throw(ME);
end

if (obj.timeIntervalType == NdgIntervalType.Constant)
    obj.timeInterval = varargin{1};
end

end