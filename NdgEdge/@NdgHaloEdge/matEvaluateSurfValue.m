function [ fM, fP ] = matEvaluateSurfValue( obj, fphys )

[ fM, fP ] = mxEvaluateSrufValue( obj.FToM, obj.FToE, obj.FToN1, obj.FToN2, fphys );
end

