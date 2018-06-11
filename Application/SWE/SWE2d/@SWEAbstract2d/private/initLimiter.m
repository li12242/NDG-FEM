function [ limiterSolver ] = initLimiter( obj )
%INITLIMITER Summary of this function goes here
%   Detailed explanation goes here

if obj.option.isKey('SWELimiterType')
    if obj.getOption('SWELimiterType') == enumSWELimiter.OnDepth
        limiterSolver = SWEDepthLimiter2d();
    elseif obj.getOption('SWELimiterType') == enumSWELimiter.OnElevation
        limiterSolver = SWEElevationLimiter2d();
    end
else
    limiterSolver = SWEDepthLimiter2d();
end
end

