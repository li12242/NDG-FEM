function matUpdateOutputResult( obj, time, fphys )

for m = 1:obj.Nmesh
    obj.outputNcFile(m).outputVar(time, [1, 2], time, fphys{m}(:,:,obj.varFieldIndex));
end

end