function setLimiterType(obj, type)
try
    obj.limiterType = NdgLimiterType( type );
catch
    msgID = 'PhysUnion:setLimiterType';
    msgtext = ['The input limiter type should be an',...
        ' enumeration from the NdgLimiterType class.'];
    ME = MException(msgID, msgtext);
    throw(ME);
end
end% func
