function [ numfluxSolver ] = initNumFluxSolver( obj )
%INITNUMFLUXSOLVER Summary of this function goes here
%   Detailed explanation goes here

if obj.option.isKey('NumFluxType')
    
    if obj.getOption('NumFluxType') == enumSWENumFlux.HLL
        numfluxSolver = SWEHLLNumFluxSolver2d( );
    elseif( obj.getOption('NumFluxType') == enumSWENumFlux.LF )
        numfluxSolver = SWELFNumFluxSolver2d( );
    elseif( obj.getOption('NumFluxType') == enumSWENumFlux.ROE )
        numfluxSolver = SWERoeNumFluxSolver2d( );
    end
else
    numfluxSolver = SWEHLLNumFluxSolver2d( );
end

end

