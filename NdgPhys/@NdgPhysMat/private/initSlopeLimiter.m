function [ limiter ] = initSlopeLimiter( physMat )
%INITSLOPELIMITER Summary of this function goes here
%   Detailed explanation goes here

if physMat.option.isKey('limiterType')
    type = physMat.getOption('limiterType');
    dim = physMat.meshUnion(1).type;
    if dim == enumMeshDim.One
        [ limiter ] = initSlopeLimiter1d( physMat.meshUnion, type );
    elseif dim == enumMeshDim.Two
        [ limiter ] = initSlopeLimiter2d( physMat.meshUnion, type );
    elseif dim == enumMeshDim.Three
        [ limiter ] = initSlopeLimiter3d( physMat.meshUnion, type );
    end
else % default non limiter
    limiter = NdgNonLimiter( physMat.meshUnion );
end

end

function [ limiter ] = initSlopeLimiter3d( mesh, type )
if ( type == enumLimiter.None )
    limiter = NdgNonLimiter( mesh );
elseif( type == enumLimiter.Vert )
    limiter = NdgVertLimiter3d( mesh );
elseif( type == enumLimiter.TVB )
    limiter = NdgTVB3d( mesh );
elseif( type == enumLimiter.BJ )
    limiter = NdgBJ3d( mesh );
end
end

function [ limiter ] = initSlopeLimiter1d( mesh, type )
if ( type == enumLimiter.None )
    limiter = NdgNonLimiter( mesh );
elseif( type == enumLimiter.Vert )
    limiter = NdgVertLimiter1d( mesh );
elseif( type == enumLimiter.TVB )
    limiter = NdgTVB1d( mesh );
elseif( type == enumLimiter.BJ )
    limiter = NdgBJ1d( mesh );
end
end

function [ limiter ] = initSlopeLimiter2d( mesh, type )
if ( type == enumLimiter.None )
    limiter = NdgNonLimiter( mesh );
elseif( type == enumLimiter.Vert )
    limiter = NdgVertLimiter2d( mesh );
elseif( type == enumLimiter.TVB )
    limiter = NdgTVB2d( mesh );
elseif( type == enumLimiter.BJ )
    limiter = NdgBJ2d( mesh );
end
end

