%======================================================================
%> @brief check whether all the properties in the class are assigned.
%>
%> @param obj object of Phys class
%======================================================================
%> This function is part of the NDGOM software. 
%> @author li12242, Tianjin University, li12242@tju.edu.cn
%======================================================================
function checkPhysProperties( obj )

% if ( isempty( obj.finalTime) )
%     msgID = 'Phys:checkPhysProperties';
%     msgtext = 'The finalTime is not assigned.';
%     throwError(msgID, msgtext);
% end

if ( isempty(obj.limiterType) ) || ( ~isa(obj.limiterType, 'NdgLimiterType') )
    msgID = 'NdgSolver:checkPhysProperties';
    msgtext = 'The limiterType is incorrect or not assigned.';
    throwError(msgID, msgtext);
end

if ( isempty(obj.temporalDiscreteType) ) || ...
        ( ~isa(obj.temporalDiscreteType, 'NdgTemporalDiscreteType') )
    msgID = 'NdgSolver:checkPhysProperties';
    msgtext = 'The temporalDiscreteType is incorrect or not assigned.';
    throwError(msgID, msgtext);
end

if ( isempty(obj.timeIntervalType) ) || ...
        ( ~isa(obj.timeIntervalType, 'NdgIntervalType') )
    msgID = 'NdgSolver:checkPhysProperties';
    msgtext = 'The timeIntervalType is incorrect or not assigned.';
    throwError(msgID, msgtext);
end

if ( obj.timeIntervalType == NdgIntervalType.Constant ) && ...
        isempty( obj.timeInterval )
    msgID = 'NdgSolver:checkPhysProperties';
    msgtext = 'The timeInterval is not assigned.';
    throwError(msgID, msgtext);
end

if ( isempty( obj.obcType ) ) || ( ~isa(obj.obcType, 'NdgBCType') )
    msgID = 'NdgSolver:checkPhysProperties';
    msgtext = 'The obcType is incorrect or not assigned.';
    throwError(msgID, msgtext);
end

if ( obj.obcType == NdgBCType.Function ) || ( obj.obcType == NdgBCType.File )
    if ( isempty( obj.obcIntervalType) ) ...
            || ( ~isa( obj.obcIntervalType, 'NdgIntervalType') )
        msgID = 'NdgSolver:checkPhysProperties';
        msgtext = 'The obcIntervalType is incorrect or not assigned.';
        throwError(msgID, msgtext);
    end
end

if ( obj.obcType == NdgBCType.File )
    if ( isempty( obj.obcFileName) )
        msgID = 'NdgSolver:checkPhysProperties';
        msgtext = 'The obcFileName is not assigned.';
        throwError(msgID, msgtext);
    end
end

fprintf('NdgSolver: All properties are well assigned.\n');
end

function throwError(msgID, msgtext)
ME = MException(msgID, msgtext);
throw(ME);
end

