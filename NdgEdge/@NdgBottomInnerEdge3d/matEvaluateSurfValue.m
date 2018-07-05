function [ fM, fP ] = matEvaluateSurfValue( obj, fphys )
%MATEVALUATESURFVALUE Summary of this function goes here
%   Detailed explanation goes here

m = obj.FToM;
[ fM, fP ] = mxEvaluateSurfValue( obj.FToE, obj.FToN1, obj.FToN2, fphys{m}  );
end

