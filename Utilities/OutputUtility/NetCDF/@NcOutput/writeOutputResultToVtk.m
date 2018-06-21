function writeOutputResultToVtk( obj, step, fieldId )
%WRITEOUTPUTRESULTTOVTK Summary of this function goes here
%   Detailed explanation goes here

field = obj.readOutputResult( step );
obj.vtkOutput.outputStepResult( step, field(:, :, fieldId) );

end

