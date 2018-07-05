function initOutput( obj, mesh2d, mesh3d )
%INITOUTPUT Summary of this function goes here
%   Detailed explanation goes here

obj.outputFile = LSWEOutput3d( obj.casename, obj.Nfield2d, obj.outputTimeInterval );
obj.outputFile.initFromMesh( mesh2d, mesh3d );

end
