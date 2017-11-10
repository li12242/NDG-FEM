function matEvaluateRHS( obj, fphys )

obj.matEvaluateCellInnerFlux( fphys );
obj.matEvaluateCellSurfFlux( fphys, obj.fext );

for m = 1:obj.Nmesh
    obj.frhs{m} = obj.frhs{m}./obj.meshUnion(m).LAVToCV;
end

end
