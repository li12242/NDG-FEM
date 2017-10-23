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
function setTemporalDiscreteType(obj, type)
% set temporalDiscreteType
try
    obj.temporalDiscreteType = NdgTemporalDiscreteType( type );
catch
    msgID = 'PhysUnion:setTemporalDiscreteType';
    msgtext = ['The input temporal discrete type should be an',...
        ' enumeration from the NdgTemporalDiscreteType class.'];
    ME = MException(msgID, msgtext);
    throw(ME);
end

end% func