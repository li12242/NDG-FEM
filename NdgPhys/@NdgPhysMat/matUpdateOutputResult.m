function matUpdateOutputResult( obj, time, step, fphys )



switch obj.getOption('outputIntervalType')
    case NdgIntervalType.Constant %> output the result by delta time
        outputResultByDeltaTime( obj, time, fphys );
        
    case NdgIntervalType.DeltaTime
        outputResultByDeltaTime( obj, time, fphys );
        
    case NdgIntervalType.DeltaStep %> output the result by delta step
        outputResultByDeltaStep( obj, step, time, fphys )
end

end

function outputResultByDeltaStep( obj, step, time, fphys )
isOutput = ( step > (obj.outputStep * obj.outputStepInterval) ) | ( abs( time - obj.ftime ) < 1e-6 );

if isOutput
    for m = 1:obj.Nmesh
        obj.outputNcFile(m).outputVar(1, time, obj.outputStep);
        obj.outputNcFile(m).outputVar(2, fphys{m}, obj.outputStep);
    end
    obj.outputStep = obj.outputStep + 1;
end
end

function outputResultByDeltaTime( obj, time, fphys )
isOutput = ( time > (obj.outputStep * obj.outputTimeInterval) ) | ( abs( time - obj.ftime ) < 1e-6 );

if isOutput
    for m = 1:obj.Nmesh
        obj.outputNcFile(m).outputVar(1, time, obj.outputStep);
        obj.outputNcFile(m).outputVar(2, fphys{m}, obj.outputStep);
    end
    obj.outputStep = obj.outputStep + 1;
end

end
