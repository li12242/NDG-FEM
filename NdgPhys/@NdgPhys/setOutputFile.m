function setOutputFile( obj, casename, intervalType, interval )

obj.outputNetcdfFileName = casename;

% set output interval type
try
    obj.outputIntervalType = NdgIntervalType( intervalType );
catch
    msgID = 'PhysUnion:setOutputFile';
    msgtext = ['The input intervalType should be an',...
        ' enumeration from the NdgIntervalType class.'];
    ME = MException(msgID, msgtext);
    throw(ME);
end

obj.outputInterval = interval;
end

