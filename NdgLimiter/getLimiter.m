function [ limiter ] = getLimiter( limiterType, meshUnion )

switch limiterType
    case NdgLimiterType.None 
        % return a limiter with doing nothing
        limiter = NdgNonLimiter(meshUnion);
        
    case NdgLimiterType.VertexLimiter 
        % return a vertex-based limiter  
        if isa( meshUnion, 'NdgMesh2d' )
            limiter = NdgVertLimiter2d( meshUnion );
        elseif isa( meshUnion, 'NdgMesh1d' )
        end
        
    otherwise
        msgID = [mfilename, ':inputLimiterTypeError'];
        msgtext = 'The input limiter type is unknown.';
        throw( MException(msgID, msgtext) );
end

end

