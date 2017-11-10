function matCloseOutputNetcdfFile( obj )
for n = 1:obj.Nmesh
    obj.outputNcFile(n).delete;
end
end